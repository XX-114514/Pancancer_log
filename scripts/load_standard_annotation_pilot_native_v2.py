"""Same-roster integer counts plus independently accepted adopted RNA readouts.

HVG ranking is left to each training fold. Candidate measurement availability is
frozen metadata, not expression from the held-out study. Counts stay in RAM.
"""
import json,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from native_archive_collection_v1 import NativeArchiveCollection
from native_annotation_review_v1 import native_measured_genes
from native_panel_aggregate_v1 import gather_csr_rows

def load_pilot(pilot_root, manifest, checkpoint, block_size, *, discovery, native_auditor):
    D=Path(discovery);A=Path(native_auditor)
    accepted=A/'outputs/native_evaluation_full_replay_acceptance_v1.json'
    while not accepted.exists():
        t=A/'provenance/supervisor_terminal_v1.json'
        if t.exists():raise ValueError('Native readout audit did not provide acceptance')
        time.sleep(5)
    acceptance=json.loads(accepted.read_text());assert acceptance['all_raw_total_and_native_RNA_readouts_replayed']
    pc=json.loads((D/'config/model_spec_v1.json').read_text());assert sha256_file(manifest)==pc['bindings']['native_manifest']['sha256']
    for item in pc['bindings'].values():assert sha256_file(item['path'])==item['sha256']
    pilot=Path(pilot_root);spec=json.loads((pilot/'config/pilot_spec_v1.json').read_text());inventory=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    selection=json.loads((D/'outputs/discovery_feature_selection_v1.json').read_text());genes=selection['candidate_gene_universe']
    contract=json.loads((D/'outputs/evaluation_RNA_contract_v1.json').read_text())
    with np.load(D/'outputs/private_shared_evaluation_vectors_v1.npz',allow_pickle=False) as f:vectors={k:f[k] for k in f.files}
    with np.load(D/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:axis={k:f[k] for k in f.files}
    reader=FullContextAtlas(manifest);native=NativeArchiveCollection(pc['native_archive_manifest']);rawpieces=[];logpieces=[];sources=[];fine=[];totals=[]
    for entry in inventory:
        ds=entry['dataset'];p=pilot/entry['selector_file'];assert sha256_file(p)==entry['selector_sha256'];positions=np.load(p,allow_pickle=False)
        measured,columns=native_measured_genes(native,ds);by=dict(zip(measured,map(int,columns)));cols=[by[g] for g in genes];assert len(cols)==len(set(cols))
        z=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];counts=[];logs=[];localtotals=[]
        for low in range(0,len(positions),block_size):
            raw=gather_csr_rows(z[key],reader.rows(ds)[positions[low:low+block_size]])
            assert np.isfinite(raw.data).all() and raw.data.min(initial=0)>=0 and np.equal(raw.data,np.floor(raw.data)).all() and raw.data.max(initial=0)<2**24
            total=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();value=raw[:,cols].astype(np.float32).tocsr();counts.append(value)
            scale=np.zeros(len(total),np.float32);ok=total>0;scale[ok]=10000/total[ok];logged=value.multiply(scale[:,None]).tocsr();logged.data=np.log1p(logged.data);logs.append(logged);localtotals.append(total)
        rawpieces.append(sp.vstack(counts,format='csr'));logpieces.append(sp.vstack(logs,format='csr'));totals.append(np.concatenate(localtotals));sources.extend([ds]*len(positions));fine.extend(reader.frozen_labels(ds)['fine'][positions].astype(str));native._fields.clear();reader._fields.clear()
        checkpoint('training_holdout_native_counts_loaded',source=ds,n_rows=len(positions))
    source=np.asarray(sources);total=np.concatenate(totals);np.testing.assert_array_equal(source,axis['source']);np.testing.assert_array_equal(total,vectors['native_total'])
    folds={r['dataset']:r['fold_group'] for r in spec['sources']}
    return dict(x=sp.vstack(logpieces,format='csr'),counts=sp.vstack(rawpieces,format='csr'),genes=genes,source=source,route=axis['inherited_route'],fine=np.asarray(fine),fold=np.asarray([folds[s] for s in source]),native_total=total,
        state={k:vectors['state_'+k] for k in contract['state']},state_genes={k:v['used_gene_names'] for k,v in contract['state'].items()},
        native_flags=vectors['native_coherent_flags'],source_assessed=contract['source_assessed'],gamma_RNA_candidate=vectors['gamma_RNA_candidate'],
        accepted_native_readout_SHA256=sha256_file(accepted),counts_origin='raw_integer_native',RNA_readout_origin='adopted_analysis_native_independently_full_replayed')
