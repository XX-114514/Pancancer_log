"""Complete 50-pathway exploratory interpretation; no mechanism selection."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    reference=Path(os.environ['PANCANCER_PUBLIC_HALLMARK_ROOT'])
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1';out=crc/'private_all50_hallmark_diagnostic_v1'
    assert not out.exists();out.mkdir()
    provenance=json.loads((reference/'reference_provenance.json').read_text());assert len(provenance['records'])==50
    # Freeze complete scope before reading learned signatures.
    (out/'plan_before_results.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      family='50 references*2 endpoints*3 directions=300, one BH family',
      full_collection_not_min_p_selection=True,interpretation_only_not_mechanism=True),indent=2)+'\n')
    sets={}
    for r in provenance['records']:
        p=reference/(r['name']+'.json');assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
        sets[r['name']]=set(json.loads(p.read_text())[r['name']]['geneSymbols'])
    signatures=pd.read_csv(crc/'discovery_active_gene_signatures.tsv',sep='\t')
    rows=[]
    for endpoint,group in signatures.groupby('cell_type'):
        bg=set(pd.read_csv(crc/'private_inputs'/(endpoint+'_discovery_genes.tsv'),header=None)[0])
        for direction in ['all','up','down']:
            genes=set(group.genes if direction=='all' else group.loc[group.direction.eq(direction),'genes'])
            assert genes.issubset(bg)
            for name,ref in sets.items():
                measured=ref&bg;overlap=genes&measured
                p=float(hypergeom.sf(len(overlap)-1,len(bg),len(measured),len(genes)))
                rows.append(dict(endpoint=endpoint,direction=direction,reference=name,background_genes=len(bg),
                  signature_genes=len(genes),reference_measurable=len(measured),overlap_genes=len(overlap),p=p))
    result=pd.DataFrame(rows);assert len(result)==300
    result['BH_all300']=multipletests(result.p,method='fdr_bh')[1]
    result.to_csv(out/'all300_private_function_statistics.tsv',sep='\t',index=False)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ALL50_FUNCTION_INTERPRETATION_COMPLETE_NOT_MECHANISM_VALIDATION',
      finished_UTC=datetime.now(timezone.utc).isoformat(),reference_sha256=hashlib.sha256((reference/'reference_provenance.json').read_bytes()).hexdigest()),indent=2)+'\n')
    print('ALL50_FUNCTION_INTERPRETATION_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
