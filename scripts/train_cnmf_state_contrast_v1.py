"""Official unmodified cNMF contrast, complete prespecified K/repeat grid.
Raw-count variance scaling is the official algorithm, distinct from Spectra.
Native-fullRNA TPM is supplied for HVG selection and gene spectra.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc,importlib.metadata
import numpy as np,anndata as ad
from cnmf import cNMF


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_CNMF_PLAN']);plan=json.loads(P.read_text())
    assert sha(__file__)==plan['producer_sha256'] and importlib.metadata.version('cnmf')=='1.7.0'
    I=R/'state_singlecell_inputs_acceptance_20261014_v2';assert sha(I/'accepted_input_inventory.json')==plan['accepted_input_inventory_sha256']
    inv=json.loads((I/'accepted_input_inventory.json').read_text());file=Path(inv['training_file']['path']);assert sha(file)==inv['training_file']['sha256']
    data=ad.read_h5ad(file);assert data.n_obs==12000 and data.obs_names.is_unique and not data.obs.gse_id.eq('GSE183904').any()
    O=R/'cNMF_state_contrast_20261014_v1';assert not O.exists();O.mkdir();summaries=[]
    for i,major in enumerate(['T cell','Macrophage']):
        name='T_cell' if major=='T cell' else 'Macrophage';root=O/name;root.mkdir();raw=data[data.obs.major.eq(major)].copy();assert raw.n_obs==6000
        raw.write_h5ad(root/'private_raw_shared_counts.h5ad',compression='gzip',compression_opts=4)
        tpm=raw.copy();tpm.X=raw.X.astype(np.float64).multiply((1000000/raw.obs.full_RNA_total.to_numpy())[:,None]).tocsr()
        tpm.write_h5ad(root/'private_native_fullRNA_TPM.h5ad',compression='gzip',compression_opts=4)
        del tpm;gc.collect()
        model=cNMF(output_dir=str(root),name='cnmf')
        model.prepare(counts_fn=str(root/'private_raw_shared_counts.h5ad'),tpm_fn=str(root/'private_native_fullRNA_TPM.h5ad'),components=plan['K'],n_iter=plan['repeats_per_K'],seed=plan['seed']+i,num_highvar_genes=plan['HVG_per_major'],max_NMF_iter=plan['max_NMF_iter'],init='random',beta_loss='frobenius')
        print(json.dumps(dict(major=major,phase='FACTORIZING_ALL_K_REPEATS',K=plan['K'],repeats=plan['repeats_per_K'],at_UTC=datetime.now(timezone.utc).isoformat())),flush=True)
        model.factorize(worker_i=0,total_workers=1,skip_completed_runs=False)
        model.combine(components=plan['K'],skip_missing_files=False)
        for k in plan['K']:
            model.consensus(k=k,density_threshold=plan['density_threshold'],local_neighborhood_size=.3,show_clustering=False)
            print(json.dumps(dict(major=major,K=k,phase='CONSENSUS_COMPLETE_PENDING_STABILITY_COMPARISON',at_UTC=datetime.now(timezone.utc).isoformat())),flush=True)
        summaries.append(dict(major=major,cells=raw.n_obs,K=plan['K'],repeats_per_K=plan['repeats_per_K'],all_requested_K_completed=True))
        (O/'progress.json').write_text(json.dumps(dict(completed_majors=len(summaries),summary=summaries),indent=2)+'\n');del raw,model;gc.collect()
    j=dict(status='PRODUCED_ALL_CNMFK_REPEATS_CONSENSUS_PENDING_COMPARISON',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),cnmf_version='1.7.0',summary=summaries,no_prior_gene_sets_used=True,normalization='Official raw-count gene-variance normalization; external native-fullRNA TPM for HVG/spectrum stages; normalized and absolute usage to be compared, not equated to Spectra amplitudes',formal_state_labels=False,Spectra_matching_complete=False,UCell_complete=False,scientific_release=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}),flush=True)


if __name__=='__main__':main()
