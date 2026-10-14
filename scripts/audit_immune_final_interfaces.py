#!/usr/bin/env python3
"""Read-only audit of released annotation, CNV, graph and CM interfaces."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from audit_immune_CM_inputs import main as recount_CM


def audit(a):
    m=a.run/'manifests'
    proof=json.loads((m/'immune_authority_release_public.json').read_text())
    cnv=json.loads((m/'CNV_logspace_v2_public.json').read_text())
    assert proof['all_checks_pass'] and proof['final_copy_all_file_checksums_verified']
    final=a.run/'outputs/immune_validation_authority_v1.zarr'
    assert not any(p.is_symlink() for p in final.rglob('*')), 'External storage links'
    g=zarr.open_group(str(final),mode='r')
    obs=read_elem(g['obs']);var=read_elem(g['var']);uns=read_elem(g['uns'])
    assert len(obs)==2670057 and len(var)==49563
    assert obs.index.is_unique and var.index.is_unique
    assert obs.dataset.nunique()==9 and obs.sampleID.nunique()==400
    assert not uns['independent_validation_frozen'] and not proof['malignancy_gold_standard']
    catalog=pd.DataFrame(json.loads(uns['subcluster_catalog_json']))
    assert len(catalog)==270 and catalog.subcluster.is_unique
    assert set(catalog.subcluster)==set(obs.subcluster.astype(str))
    actual=obs.subcluster.astype(str).value_counts()
    size_column=next(k for k in ['n_cells','n_atlas_cells'] if k in catalog)
    assert all(actual[r.subcluster]==getattr(r,size_column) for r in catalog.itertuples())
    assert obs.groupby('subcluster',observed=True).major.nunique().eq(1).all()
    statuses=obs.malignant_status.astype(str).value_counts().to_dict()
    assert statuses==cnv['cell_malignant_status_counts']==proof['malignancy_status_counts']
    assert cnv['n_method_failed']==0 and cnv['all_sample_IDs_exact'] and cnv['all_cell_IDs_exact']
    assert cnv['old_raw_integer_CNVs_scientifically_quarantined']
    yes=lambda x:x.astype(str).str.lower().isin(['true','1'])
    evaluated=yes(obs.cnv_evaluable)
    assert not obs.loc[evaluated,'malignant_status'].astype(str).str.startswith('not_evaluated').any()
    assert obs.loc[obs.malignant_status.astype(str).eq('CNV_supported_malignant'),'cnv_evaluable'].astype(str).str.lower().isin(['true','1']).all()
    total=pd.to_numeric(obs.CNV_native_total_counts,errors='coerce')
    assert np.isfinite(total[evaluated]).all() and total[evaluated].gt(0).all()
    graphs={}
    for name in g['obsp']:
        graph=g['obsp'][name]
        assert list(graph.attrs['shape'])==[len(obs),len(obs)]
        assert len(graph['indptr'])==len(obs)+1
        previous=0
        for lo in range(0,len(graph['indptr']),1000000):
            ptr=np.asarray(graph['indptr'][lo:lo+1000000])
            assert ptr[0]>=previous and np.all(np.diff(ptr)>=0)
            previous=int(ptr[-1])
        assert previous==len(graph['data'])==len(graph['indices'])
        for lo in range(0,previous,1000000):
            values=np.asarray(graph['data'][lo:lo+1000000])
            indices=np.asarray(graph['indices'][lo:lo+1000000])
            assert np.isfinite(values).all() and (values>=0).all()
            assert (indices>=0).all() and (indices<len(obs)).all()
        graphs[name]={'n_rows':len(obs),'n_stored_edges':previous}
    assert 'connectivities' in graphs and 'distances' in graphs
    cm_reports=json.loads(uns['CM_scope_reports_json'])
    cm_memberships=json.loads(uns['CM_node_membership_json'])
    cm_index=json.loads((m/'CM_recovery_v3_public.json').read_text())
    assert len(cm_reports)==15 and set(cm_memberships)=={r['scope'] for r in cm_index['results']}
    for entry in cm_reports:
        discovery=entry['discovery'];name=entry['input_scope']['scope']
        accepted=next(r for r in cm_index['results'] if r['scope']==name)
        recorded=json.loads((a.run/'CM'/name/accepted['discovery_attempt']/'CM_public.json').read_text())
        assert discovery==recorded
        member_path=a.run/'CM'/name/accepted['discovery_attempt']/'each_node.tsv'
        try:
            member_frame=pd.read_csv(member_path,sep='\t')
            members=member_frame.astype(object).where(pd.notna(member_frame),None).to_dict(orient='records')
        except pd.errors.EmptyDataError:
            assert discovery['n_network_nodes']==0
            members=[]
        assert cm_memberships[name]==members
    recount_CM(argparse.Namespace(run=a.run,scratch=a.scratch,observation_zarr=final,
        report=m/'CM_final_archive_recount_v1_public.json'))
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),status='PASS_FINAL_ARCHIVE_INTERFACES',
        n_cells=len(obs),n_features=len(var),n_sources=9,n_original_libraries=400,n_subclusters=270,
        standalone_final_archive=True,annotation_catalog_exact=True,CNV_exact_mapping_and_statuses=True,
        CNV_native_total_counts_storage_dtype=str(obs.CNV_native_total_counts.dtype),graphs=graphs,
        n_CM_scopes=15,CM_discovery_reports_exact=True,CM_input_frequencies_recounted_from_final_archive=True,
        malignancy_status_counts=statuses,independent_clinical_validation=False,malignancy_gold_standard=False)
    with (m/'immune_final_interface_audit_v1_public.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--scratch',type=Path,required=True)
    audit(p.parse_args())
