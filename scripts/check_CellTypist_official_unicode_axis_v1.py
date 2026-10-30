"""Compare ordinary and Unicode-metadata official predictions, not accuracy."""
import argparse,copy,datetime,hashlib,inspect,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
    cfg=load(R/'config/spec_v1.json');old=Path(cfg['existing_classifier']);oc=load(old/'config/spec_v1.json');H=Path(oc['model']);J=Path(oc['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(old/'scripts'));sys.path.insert(0,str(R/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes,model_fingerprint
    from native_panel_aggregate_v1 import gather_csr_rows
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    from CellTypist_official_unicode_axis_v1 import official_annotate_unicode_axis,independent_exact_membership_check
    import anndata,celltypist,logging
    from celltypist.models import Model
    from celltypist.classifier import Classifier
    logging.getLogger('celltypist').setLevel(logging.WARNING);storage_gate(R,cfg);assert celltypist.__version__=='1.7.1'
    source_hash=hashlib.sha256(inspect.getsource(Classifier.celltype).encode()).hexdigest();assert source_hash==cfg['official_celltype_source_SHA256']
    # Metadata counterexamples cover prefix collisions, long placeholders and
    # Unicode without truncation or alias matching.
    q=['A','AA','A-1','__UNRESOLVED_NATIVE_COLUMN_'+('0'*120),'\u03b3_gene'];f=['AA','\u03b3_gene'];np.testing.assert_array_equal(independent_exact_membership_check(q,f),[False,True,False,False,True])
    with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as z:axes={k:z[k] for k in z.files}
    lineage=load(H/'outputs/source_input_lineage_v1.json');native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);probes=[]
    for source in cfg['probe_sources']:
        code=next(i for i,r in enumerate(lineage) if r['source']==source);r=lineage[code];ix=np.flatnonzero(axes['source_code']==code)[:cfg['n_probe_cells']];store=native.store(source);genes,columns=native_measured_genes(native,source);by=dict(zip(columns.tolist(),genes));names=[by.get(j,'__UNRESOLVED_NATIVE_COLUMN_'+str(j)) for j in range(int(store[r['counts_key']].attrs['shape'][1]))]
        x=gather_csr_rows(store[r['counts_key']],axes['native_row'][ix],max_span=8192).astype(np.float64).tocsr()
        if 'ambient_removed_counts' in store['layers']:x=(x-gather_nonnegative_float_delta(store['layers/ambient_removed_counts'],axes['native_row'][ix])).tocsr()
        x.eliminate_zeros();assert x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();tot=np.asarray(x.sum(axis=1)).ravel();assert (tot>0).all();x=x.multiply((10000/tot)[:,None]).tocsr();x.data=np.log1p(x.data);ad=anndata.AnnData(x,var=pd.DataFrame(index=names))
        for key in cfg['probe_models']:
            assert sha(oc['models'][key]['path'])==oc['models'][key]['sha256'];model=Model.load(oc['models'][key]['path']);before=model_fingerprint(model);mask=independent_exact_membership_check(names,model.features)
            if not mask.any():continue
            start=time.monotonic();ordinary=celltypist.annotate(ad.copy(),model=copy.deepcopy(model),mode='best match',majority_voting=False);slow=time.monotonic()-start
            start=time.monotonic();fast=official_annotate_unicode_axis(ad.copy(),model);elapsed=time.monotonic()-start
            for field in ['decision_matrix','probability_matrix']:
                a=getattr(ordinary,field);b=getattr(fast,field);assert a.index.equals(b.index) and a.columns.equals(b.columns);np.testing.assert_array_equal(a.to_numpy(),b.to_numpy())
            np.testing.assert_array_equal(ordinary.predicted_labels.predicted_labels.astype(str).to_numpy(),fast.predicted_labels.predicted_labels.astype(str).to_numpy());assert model_fingerprint(model)==before
            probes.append(dict(source=source,model=key,n_probe_cells=len(ix),n_exact_matching_features=int(mask.sum()),decision_and_probability_bitwise_identical=True,labels_identical=True,original_model_unchanged=True,ordinary_seconds=slow,unicode_metadata_seconds=elapsed,biological_scope_or_accuracy_not_tested=True))
        native._fields.clear()
    assert len(probes)==len(cfg['probe_sources'])*len(cfg['probe_models']);assert hashlib.sha256(inspect.getsource(Classifier.celltype).encode()).hexdigest()==source_hash
    accept=dict(schema='accepted_official_CellTypist_unicode_gene_metadata_equivalence_v1',celltypist_version=celltypist.__version__,official_celltype_source_SHA256=source_hash,helper_SHA256=sha(R/'scripts/CellTypist_official_unicode_axis_v1.py'),probes=probes,expression_model_numbers_and_official_functions_unchanged=True,numerical_and_gene_membership_counterexamples_pass=True,no_global_numpy_or_package_monkeypatch=True,majority_voting_or_scope_not_changed=True,not_a_biological_annotation_accuracy_test=True)
    save(R/'outputs/official_unicode_axis_equivalence_acceptance_v1.json',accept);storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='OFFICIAL_CELLTYPIST_UNICODE_GENE_METADATA_EQUIVALENCE_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/official_unicode_axis_equivalence_acceptance_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
