"""Reproducible replay of public V7 aggregates; never assigns cell labels."""
from pathlib import Path
import json,hashlib,argparse
import pandas as pd
from evidence_engine import evaluate

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    rpath=a.repo/'resources/annotation/v7_panel_repair_v1/registry.json'
    r=json.loads(rpath.read_text());base=a.repo/'inventories/analysis/2026-11-05'
    old=base/'v7_panel_semantics_audit_v1'
    ep=old/'sources/expanded_panel_gene_evidence_v1.tsv.gz';gp=old/'inputs/legacy_expanded_cluster_gates.tsv'
    e=pd.read_csv(ep,sep='\t');g=pd.read_csv(gp,sep='\t').set_index('cluster_uid')
    keys=['cluster_uid','gene'];value=['n_eligible_cells','pct_in']
    assert (e.groupby(keys)[value].nunique(dropna=False)<=1).all().all(), 'Conflicting old measurements'
    e=e.drop_duplicates(keys).copy();e['n_cells']=e.cluster_uid.map(g.n_cells)
    ex=e.rename(columns={'cluster_uid':'cluster','n_eligible_cells':'n_measured','pct_in':'pct_expressing_measured'})[['cluster','gene','n_cells','n_measured','pct_expressing_measured']]
    t6=base/'v7_literature_guided_v1/T6_diagnostic_v1';tp=t6/'T6_measured_marker_expression.tsv';dp=t6/'T6_top10_RNA_DEG.tsv'
    inputs={'old40':(ex,None),'T6_12':(pd.read_csv(tp,sep='\t',dtype={'cluster':str}),pd.read_csv(dp,sep='\t',dtype={'cluster':str}))}
    results={}
    for name,(expr,de) in inputs.items():
        d=evaluate(r,expr,de=de,de_complete=False)
        assert not d.final_label_assigned.any()
        assert len(d)==expr.cluster.nunique()*len(r['panels'])
        d.replace('',pd.NA).to_csv(a.out/(name+'_all_panel_evidence.tsv.gz'),sep='\t',index=False,na_rep='NA',compression={'method':'gzip','mtime':0})
        counts=expr.groupby('cluster').n_cells.first()
        results[name]={'clusters':len(counts),'cells_in_scope':int(counts.sum()),'evidence_rows':len(d),'measurement_rows':len(expr),'DEG':'not_supplied_legacy_statistic_not_reused' if de is None else 'top10_truncated_no_absence_inference','final_label_changes':0}
    # A compact observed contrast, not a new biological label assignment.
    t=inputs['T6_12'][0];genes=['FOXP3','IL2RA','IKZF2','CD4','CD8A','CD8B','CXCL13','CXCR5','BCL6','PDCD1','GZMK']
    t[t.gene.isin(genes)].to_csv(a.out/'T6_selected_observed_markers.tsv',sep='\t',index=False,na_rep='NA')
    sources=[rpath,ep,gp,tp,dp]
    manifest={'method':'all candidates, distinct axes, no ranking or automatic relabeling','results':results,
      'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('evidence_engine.py')]},
      'inputs':{str(p.relative_to(a.repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
      'limitations':['Old40 and T6 scopes overlap: never sum their cell counts','New registry genes absent from old inputs are explicitly unmeasured/not supplied','These are aggregate replays, not new full-matrix biological fits','Top10 DEG is truncated, absence from it is not evidence of no enrichment']}
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(results))
if __name__=='__main__':main()
