"""Descriptive calibration of source-common measured score against full CRC.
No signatures, thresholds, or quality coefficients are refitted.
"""
import gzip,hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem
from scipy.io import mmread
from scipy.stats import spearmanr


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    out=crc/'private_observed_part_calibration_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_results.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      source_common_axes='all4 nonCRC native sources',scope='both endpoints discovery and holdout',
      no_refitting=True,no_new_significance_tests=True,interpretation='calibration not independent biological validation'),indent=2)+'\n')
    axes=[]
    for source in ['GSE148071','GSE161529','GSE182109','GSE183904']:
        with h5py.File(root/'state_singlecell_inputs_20261014_v1/private_native_sources'/(source+'.h5ad')) as f:
            axes.append(set(read_elem(f['var']).index))
    common=set.intersection(*axes);rows=[]
    for endpoint in ['Tnode9','Macnode18']:
        w=pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_genes.tsv'),sep='\t')
        observed=w.gene.isin(common).to_numpy();fraction=float(np.abs(w.weight[observed]).sum()/np.abs(w.weight).sum())
        qc=pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_quality.tsv'),sep='\t').set_index('term').value
        for partition in ['discovery','holdout']:
            prefix=crc/'private_inputs'/(endpoint+'_'+partition)
            genes=pd.read_csv(str(prefix)+'_genes.tsv',header=None)[0]
            indices=pd.Index(genes).get_indexer(w.gene);assert (indices>=0).all()
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as f:mat=mmread(f).tocsr()[indices].toarray().T
            md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t')
            z=(mat-w['mean'].to_numpy())/w.sd.to_numpy()
            correction=qc['(Intercept)']+md.cellQ.to_numpy()*qc.cellQ+md.geneQ.to_numpy()*qc.geneQ
            full=z@w.weight.to_numpy()-correction
            partial=z[:,observed]@w.weight.to_numpy()[observed]-correction
            patient=pd.DataFrame({'sample':md['sample'],'full':full,'partial':partial}).groupby('sample').mean()
            rows.append(dict(endpoint=endpoint,partition=partition,absolute_weight_retained=fraction,
                             cells=len(md),patient_blocks=len(patient),cell_score_rho=float(spearmanr(full,partial).statistic),
                             patient_score_rho=float(spearmanr(patient.full,patient.partial).statistic),
                             centered_component_RMSE=float(np.sqrt(np.mean((partial-full)**2)))))
    pd.DataFrame(rows).to_csv(out/'all4_private_score_calibration.tsv',sep='\t',index=False)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='OBSERVED_PART_CALIBRATION_COMPLETE_NO_REFITTING',
      finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print('OBSERVED_PART_CALIBRATION_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
