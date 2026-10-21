"""Author subtype counts and exact cached RNA QC/quality; case data RAM only.

This is a new source-specific measurement; no marker-parent equivalence,
formal project identity update or clinical association is implied.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import resource
import numpy as np
import pandas as pd

CLINICAL=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(16*1024**2),b''):h.update(block)
    return h.hexdigest()


def aggregate(mask,case,total,genes,denominator,ncase):
    counts=np.bincount(case[mask],minlength=ncase)
    means=np.column_stack([np.bincount(case[mask],weights=np.log1p(x[mask]),minlength=ncase) for x in [total,genes]])
    means/=np.maximum(1,counts)[:,None]
    table=pd.DataFrame(dict(case=case[mask],RNA=np.log1p(total[mask]),genes=np.log1p(genes[mask])))
    independent=table.groupby('case').mean().reindex(range(ncase))
    np.testing.assert_allclose(means[counts>0],independent.to_numpy()[counts>0],atol=1e-10,rtol=0)
    quality=[[float(means[i,0]) if counts[i] else None,float(means[i,1]) if counts[i] else None,
              float(np.log1p(counts[i])),float(counts[i]/denominator[i]) if denominator[i] else None] for i in range(ncase)]
    return counts,quality


def read_source(root,plan):
    available=int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:')))*1024
    assert available>=4*1024**3,'AUTHOR_COUNT_SOURCE_MEMORY_HEADROOM_NOT_ACCEPTED'
    root=Path(root);native=root/'NSCLC_full_native_gene_units_20261016_v1'
    assert sha(native/'private_terminal.json')==plan['native_terminal_SHA256']
    terminal=json.loads((native/'private_terminal.json').read_text())
    names=['private_all_native_units_and_quality.tsv','private_cell_unit_index.npy','private_full_RNA.npy','private_detected_genes.npy']
    for name in names:assert sha(native/name)==terminal['outputs'][name]['sha256']
    units=pd.read_csv(native/names[0],sep='\t',keep_default_na=False).sort_values('unit_index')
    assert units.unit_index.tolist()==list(range(len(units)))
    assert not units.duplicated(['sampleID','sub_cell_type']).any()
    unit=np.load(native/names[1],allow_pickle=False);total=np.load(native/names[2],allow_pickle=False)
    detected=np.load(native/names[3],allow_pickle=False)
    assert len(unit)==len(total)==len(detected)==1254749 and unit.min()>=0 and unit.max()<len(units)
    np.testing.assert_array_equal(np.bincount(unit,minlength=len(units)),units.n_cells.to_numpy(int))
    metadata=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(metadata)==plan['native_metadata_SHA256']
    md=pd.read_csv(metadata,usecols=['sampleID','sub_cell_type','major_cell_type','total_counts','total_counts_mt'],
                   dtype={'sampleID':str,'sub_cell_type':str,'major_cell_type':str},keep_default_na=False)
    assert np.array_equal(md.sampleID,units.sampleID.to_numpy()[unit])
    assert np.array_equal(md.sub_cell_type,units.sub_cell_type.to_numpy()[unit])
    assert np.array_equal(md.major_cell_type,units.major_cell_type.to_numpy()[unit])
    np.testing.assert_array_equal(pd.to_numeric(md.total_counts),total)
    mito=pd.to_numeric(md.total_counts_mt).to_numpy();assert np.isfinite(mito).all() and (mito>=0).all() and (mito<=total).all()
    case_names=sorted(set(units.sampleID));assert len(case_names)==243
    case=pd.Index(case_names).get_indexer(md.sampleID);assert (case>=0).all()
    clinical_path=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv'
    assert sha(clinical_path)==plan['clinical_table_SHA256']
    clinical=pd.read_csv(clinical_path,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID')
    assert clinical.index.is_unique and set(clinical.index)==set(case_names);clinical=clinical.loc[case_names]
    numeric=clinical[CLINICAL+['response']].apply(pd.to_numeric,errors='coerce')
    eligible=clinical.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE').to_numpy() & np.isfinite(numeric.to_numpy()).all(axis=1)
    assert np.array_equal(eligible,clinical.prospective_RNA_clinical_eligible.map({'True':True,'False':False}))
    definitions=plan['author_lineage_names'];assert set(definitions)=={'Treg','CD8T','CD4T_nonTreg'}
    labels=units.sub_cell_type.to_numpy();masks={k:np.isin(labels,v)[unit] for k,v in definitions.items()}
    assert all(m.any() for m in masks.values()) and np.sum(list(masks.values()),axis=0).max()==1
    assert set(md.major_cell_type[np.any(list(masks.values()),axis=0)])=={'T/NK cell'}
    state_names=definitions['Treg']+definitions['CD8T'];assert len(state_names)==13 and len(set(state_names))==13
    for parent in definitions:
        assert set(definitions[parent])<=set(labels)
    output=[];coverage=[];source_QC_parent_totals={}
    for qc in ['native_QC500_250_MT20','all_author_retained']:
        good=(total>=500)&(detected>=250)&(total>0)&(mito<=.2*total) if qc=='native_QC500_250_MT20' else np.ones(len(unit),bool)
        chosen={k:m&good for k,m in masks.items()}
        denominator=np.bincount(case[good & np.any(list(masks.values()),axis=0)],minlength=len(case_names))
        parent_counts={};quality={}
        for k,mask in chosen.items():parent_counts[k],quality[k]=aggregate(mask,case,total,detected,denominator,len(case_names))
        source_QC_parent_totals[qc]={k:int(v.sum()) for k,v in parent_counts.items()}
        state_counts=[]
        for label in state_names:
            mask=good & (labels[unit]==label)
            counts=np.bincount(case[mask],minlength=len(case_names))
            independent=pd.Series(case[mask]).value_counts().reindex(range(len(case_names)),fill_value=0).to_numpy(int)
            np.testing.assert_array_equal(counts,independent);state_counts.append(counts)
        state_counts=np.asarray(state_counts).T
        np.testing.assert_array_equal(state_counts[:,:3].sum(axis=1),parent_counts['Treg'])
        np.testing.assert_array_equal(state_counts[:,3:].sum(axis=1),parent_counts['CD8T'])
        for i in range(len(case_names)):
            output.append(dict(case_index=i,QC=qc,author_counts=state_counts[i].tolist(),
                parent_counts={k:int(v[i]) for k,v in parent_counts.items()},quality={k:v[i] for k,v in quality.items()},
                clinical=[float(v) if np.isfinite(v) else None for v in numeric[CLINICAL].iloc[i]],
                response=float(numeric.response.iloc[i]) if np.isfinite(numeric.response.iloc[i]) else None,
                clinical_eligible=bool(eligible[i]),cancer_type=clinical.cancer_type.iloc[i],drug=clinical['anti-PD1_therapy'].iloc[i]))
        for cancer in ['ALL','LUAD','LUSC']:
            base=eligible if cancer=='ALL' else eligible & clinical.cancer_type.eq(cancer).to_numpy()
            both=base & (parent_counts['Treg']>=20)&(parent_counts['CD8T']>=20)
            summary=dict(QC=qc,cancer=cancer,clinical_metadata_units=int(base.sum()),both_parent_ge20=int(both.sum()),
                MPR_pCR=int((both & numeric.response.eq(1).to_numpy()).sum()),nonMPR=int((both & numeric.response.eq(0).to_numpy()).sum()),
                parent_cell_totals_clinical_eligible={k:int(v[base].sum()) for k,v in parent_counts.items()},
                all13_author_states_variable={state:bool(len(np.unique(state_counts[both,j]/parent_counts['Treg' if j<3 else 'CD8T'][both]))>1) for j,state in enumerate(state_names)})
            coverage.append(summary)
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    assert peak<=2*1024**3,'AUTHOR_COUNT_SOURCE_RSS_CAP'
    return dict(rows=output,state_names=state_names,source_coverage=coverage,
        proof=dict(all_original_cell_unit_case_clinical_keys_and_SHA_verified=True,
            original_QC_exactly_applied_to_all_author_lineages=True,source_annotation_counts_not_marker_positive_relabelling=True,
            all_state_counts_and_parent_mean_log_RNA_gene_quality_independently_checked=True,
            fraction_denominator='Sum of these explicitly declared author CD4T and CD8T families, not all pan-T or all TME',
            all243_source_QC_parent_cell_totals=source_QC_parent_totals,
            original_CM_or_formal_project_identity_not_changed=True,clinical_or_coordination_effects_not_computed=True,
            source_peak_RSS_bytes=int(peak)))
