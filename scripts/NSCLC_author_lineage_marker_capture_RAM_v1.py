"""Exact single-gene detection in author lineages from accepted unit caches.

This is permitted only if QC removes no author-parent cells. Marginal counts
cannot reconstruct Boolean multi-gene coexpression or individual cell labels.
"""
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess


def read_source(root,plan):
    import numpy as np
    import pandas as pd
    root=Path(root);native=root/'NSCLC_full_native_gene_units_20261016_v1'
    def sha(path):
        h=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(16*1024**2),b''):h.update(block)
        return h.hexdigest()
    assert sha(native/'private_terminal.json')==plan['native_terminal_SHA256']
    terminal=json.loads((native/'private_terminal.json').read_text())
    registry=native/'private_all_native_units_and_quality.tsv';positive_path=native/'positive_cell_counts.npy'
    for path in [registry,positive_path]:assert sha(path)==terminal['outputs'][path.name]['sha256']
    units=pd.read_csv(registry,sep='\t',keep_default_na=False).sort_values('unit_index')
    assert units.unit_index.tolist()==list(range(len(units))) and not units.duplicated(['sampleID','sub_cell_type']).any()
    gene_path=root/'GSE243013_genes.csv.gz';assert sha(gene_path)==plan['native_gene_axis_SHA256']
    genes=pd.read_csv(gene_path,keep_default_na=False).geneSymbol.str.strip().str.upper()
    assert genes.is_unique;columns=pd.Index(genes).get_indexer(plan['markers']);assert (columns>=0).all()
    positive=np.load(positive_path,mmap_mode='r',allow_pickle=False);assert positive.shape==(len(units),len(genes))
    case_names=sorted(set(units.sampleID));assert len(case_names)==243
    cases=pd.Index(case_names).get_indexer(units.sampleID);assert (cases>=0).all()
    parent_counts={};marker_counts={}
    for parent,names in plan['author_lineage_names'].items():
        mask=units.sub_cell_type.isin(names).to_numpy()
        parent_counts[parent]=np.bincount(cases[mask],weights=units.n_cells.to_numpy()[mask],minlength=len(case_names)).astype(np.int64)
        if parent=='CD8T':continue
        values=[]
        for column in columns:
            counts=np.asarray(positive[:,column],np.int64)
            assert (counts>=0).all() and (counts<=units.n_cells.to_numpy()).all()
            total=np.bincount(cases[mask],weights=counts[mask],minlength=len(case_names)).astype(np.int64)
            independent=pd.DataFrame(dict(case=cases[mask],count=counts[mask])).groupby('case')['count'].sum().reindex(range(len(case_names)),fill_value=0).to_numpy(np.int64)
            np.testing.assert_array_equal(total,independent);assert (total<=parent_counts[parent]).all()
            values.append(total)
        marker_counts[parent]=np.asarray(values).T
    rows=[dict(case_index=i,parent_counts={k:int(v[i]) for k,v in parent_counts.items()},
          marker_counts={k:v[i].tolist() for k,v in marker_counts.items()}) for i in range(len(case_names))]
    return dict(rows=rows,proof=dict(source='GSE243013',markers=plan['markers'],
        current_positive_unit_cache_SHA256=terminal['outputs'][positive_path.name]['sha256'],
        source_parent_marker_detection_totals={k:v.sum(axis=0).tolist() for k,v in marker_counts.items()},
        exact_gene_marginal_counts_not_multigene_coexpression=True,all_case_counts_independent_groupby_checked=True))


def augment(payload,task,plan):
    proof=json.loads((task/'measurement_v1/private_measurement_proof.json').read_text())
    totals=proof['all243_source_QC_parent_cell_totals']
    assert totals['native_QC500_250_MT20']==totals['all_author_retained'], 'QC_IDENTITY_REQUIRED_FOR_EXACT_MARGINAL_REUSE'
    source_plan=plan['marker_source_plan'];code=Path(__file__).read_text()
    program=code+'\ntry:\n result=read_source('+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+','+repr(source_plan)+')\n print(json.dumps(dict(kind="RAM_payload",value=result),allow_nan=False),flush=True)\n'
    program+='except Exception as error:\n print(json.dumps(dict(kind="failed",error_type=type(error).__name__)),flush=True)\n raise SystemExit(1)\n'
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15']
    if os.environ.get('PANCANCER_READONLY_PROXY_COMMAND'):args+=['-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND']]
    args += [os.environ['PANCANCER_READONLY_SSH_TARGET'],
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -']
    result=subprocess.run(args,input=program,capture_output=True,text=True,timeout=120)
    assert result.returncode==0,'EXACT_SOURCE_MARKER_CACHE_NOT_ACCEPTED'
    message=json.loads(result.stdout);assert message['kind']=='RAM_payload'
    markers=message['value'];lookup={r['case_index']:r for r in markers['rows']};assert len(lookup)==243
    for row in payload['rows']:
        ref=lookup[row['case_index']]
        assert row['parent_counts']==ref['parent_counts']
        row['marker_counts']=ref['marker_counts']
    folder=task/'marker_measurement_v1';folder.mkdir(exist_ok=False)
    markers['proof'].update(QC_parent_set_identity_exactly_qualified=True,all_RAM_case_parent_count_alignment_checked=True,
        case_marker_counts_RAM_only=True,formal_project_identity_or_old_family_not_changed=True)
    (folder/'private_source_proof.json').write_text(json.dumps(markers['proof'],indent=2)+'\n')
    return payload
