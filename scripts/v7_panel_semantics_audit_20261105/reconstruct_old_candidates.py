"""Read-only audit of historical aggregate evidence; does not assign cells."""
from pathlib import Path
import ast, json, hashlib, runpy
import numpy as np
import pandas as pd

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-dir',type=Path,required=True)
parser.add_argument('--output-dir',type=Path,required=True)
args=parser.parse_args()
HERE=args.input_dir
PRIOR=HERE/'inputs'
OUT=args.output_dir
OUT.mkdir(parents=True,exist_ok=False)
spec_module=runpy.run_path(str(HERE/'sources/pancancer_evidence_marker_specs.py'))
specs=spec_module['SPECS']; blacklisted=spec_module['is_blacklisted_marker']
var=set(json.loads((PRIOR/'legacy_panel_var.json').read_text()))
aliases={'KIM1':'HAVCR1'}
catalog=[]
for major,s in specs.items():
    for name,genes in s.panels.items():
        effective=[aliases.get(g,g) for g in genes]
        canon=[aliases.get(g,g) for g in s.canonical[name]]
        catalog.append(dict(major=major,old_panel=name,genes=genes,canonical=s.canonical[name],
            missing_from_legacy_object=[g for g in effective if g not in var],
            blacklisted_genes=[g for g in effective if blacklisted(g)],
            canonical_not_in_panel=sorted(set(canon)-set(effective)),
            min_pct=s.min_pct_in,min_log_difference=s.min_log2fc,min_delta_pct=s.min_delta_pct,
            nominal_resolution=s.resolution))
(OUT/'panel_catalogue.json').write_text(json.dumps(catalog,indent=2)+'\n')
major_source=ast.parse((HERE/'sources/02_major_evidence_v3.py').read_text())
major={n.target.id:ast.literal_eval(n.value) for n in major_source.body if isinstance(n,ast.AnnAssign) and isinstance(n.target,ast.Name) and n.target.id in ['MAJOR_MARKERS','MAJOR_EXCLUDE_MARKERS']}
(OUT/'major_panel_catalogue.json').write_text(json.dumps(major,indent=2)+'\n')
gene=pd.read_csv(HERE/'sources/expanded_panel_gene_evidence_v1.tsv.gz',sep='\t')
old=pd.read_csv(PRIOR/'legacy_expanded_cluster_gates.tsv',sep='\t')
initial=pd.read_csv(PRIOR/'legacy_initial_cluster_gates.tsv',sep='\t')
cross=pd.read_csv(PRIOR/'v7_label_legacy_cluster_crosswalk.tsv',sep='\t')
assert not gene.duplicated(['cluster_uid','subtype','gene']).any()
calc=gene.pct_in*np.maximum(gene.log2fc,0)*np.maximum(gene.delta_pct,0)
assert np.allclose(calc,gene.gene_evidence_score,atol=1e-9)
assert np.allclose(gene.log2fc,(gene.mean_in-gene.mean_out)/np.log(2),atol=1e-9)
for maj,g in gene.groupby('major_celltype'):
    s=specs[maj]
    passes=(g.pct_in>=s.min_pct_in)&(g.log2fc>=s.min_log2fc)&(g.delta_pct>=s.min_delta_pct)&(g.n_eligible_cells>0)
    assert (passes==g.passes_threshold).all()
rows=[]
for _,o in old.iterrows():
    gs=gene[gene.cluster_uid==o.cluster_uid]
    candidates=[]
    for name in specs[o.major_celltype].panels:
        q=gs[gs.subtype==name]; passed=q[q.passes_threshold]; cp=passed[passed.is_canonical]
        score=float(passed.gene_evidence_score.mean()) if len(passed) else -np.inf
        candidates.append(dict(cluster_uid=o.cluster_uid,major=o.major_celltype,panel=name,score=score,
            n_pass=len(passed),n_canonical_pass=len(cp),pass_genes=','.join(passed.gene),canonical_pass_genes=','.join(cp.gene),
            minimal_gene_gate=len(passed)>=2 and len(cp)>=1,
            old_contamination_name='contamination' in name.lower(),
            v7_n_cells=int(cross.loc[cross.cluster_uid==o.cluster_uid,'n_cells'].sum()),
            old_rescued=bool(o.rescued)))
    candidates.sort(key=lambda x:x['score'],reverse=True)
    assert candidates[0]['panel']==o.top_subtype
    expected=o.top_score
    assert (not np.isfinite(candidates[0]['score']) and not np.isfinite(expected)) or np.isclose(candidates[0]['score'],expected,atol=1e-8)
    second=candidates[1]
    assert second['panel']==o.second_subtype
    assert (not np.isfinite(second['score']) and not np.isfinite(o.second_score)) or np.isclose(second['score'],o.second_score,atol=1e-8)
    margin=(candidates[0]['score']-(second['score'] if np.isfinite(second['score']) else 0.0)) if np.isfinite(candidates[0]['score']) else np.nan
    assert (pd.isna(margin) and pd.isna(o.margin)) or np.isclose(margin,o.margin,atol=1e-8)
    assert candidates[0]['n_pass']==o.n_full_panel_pass
    assert candidates[0]['n_canonical_pass']==o.n_canonical_pass
    for i,c in enumerate(candidates,1):c['old_rank']=i
    rows+=candidates
candidates=pd.DataFrame(rows)
effects=[]
for _,o in old.iterrows():
    cs=candidates[candidates.cluster_uid==o.cluster_uid]; first=cs.iloc[0]
    alt=cs[(cs.old_rank>1)&cs.minimal_gene_gate&~cs.old_contamination_name]
    effects.append(dict(cluster_uid=o.cluster_uid,major=o.major_celltype,n_cells=int(o.n_cells),
        v7_n_cells=int(first.v7_n_cells),old_top=o.top_subtype,old_second=o.second_subtype,
        old_margin=o.margin,top_pass_genes=first.pass_genes,top_canonical_genes=first.canonical_pass_genes,
        top_fails_minimal_gene_gate=not bool(first.minimal_gene_gate),
        lower_candidates_passing_minimal_gene_gate=';'.join(alt.panel),
        lower_candidates_n=len(alt),top_has_no_passing_genes=int(first.n_pass)==0,
        contamination_top=bool(first.old_contamination_name),old_rescued=bool(o.rescued),old_reason=o.reason))
effects=pd.DataFrame(effects)
shadow=effects[effects.top_fails_minimal_gene_gate&(effects.lower_candidates_n>0)]
for name,df in [('candidate_reconstruction',candidates),('cluster_rule_audit',effects),('shadowed_gene_supported_candidates',shadow)]:
    df.to_csv(OUT/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
overlap=[]
for maj,s in specs.items():
    names=list(s.panels)
    for i,a in enumerate(names):
        for b in names[i+1:]:
            g=set(s.panels[a])&set(s.panels[b]); c=set(s.canonical[a])&set(s.canonical[b])
            if g:overlap.append(dict(major=maj,panel_a=a,panel_b=b,shared_genes=','.join(sorted(g)),shared_canonical=','.join(sorted(c))))
pd.DataFrame(overlap).to_csv(OUT/'within_lineage_panel_overlaps.tsv',sep='\t',index=False,na_rep='NA')
summary=dict(n_subtype_panels=len(catalog),n_major_positive=len(major['MAJOR_MARKERS']),n_major_exclusion=len(major['MAJOR_EXCLUDE_MARKERS']),
    n_contamination_panels=sum('contamination' in c['old_panel'] for c in catalog),
    n_ifn_cycling_panels=sum(c['old_panel'].startswith(('IFN_','Cycling_')) for c in catalog),
    n_gene_rows=len(gene),n_reconstructed_candidates=len(candidates),n_reconstructed_clusters=len(old),
    all_scores_and_top_candidates_reproduced=True,shadow_clusters=shadow.cluster_uid.tolist(),
    shadow_cells=int(shadow.v7_n_cells.sum()),no_signal_clusters=effects.loc[effects.top_has_no_passing_genes,'cluster_uid'].tolist(),
    contamination_top_clusters=effects.loc[effects.contamination_top,'cluster_uid'].tolist(),
    missing_genes={c['major']+'/'+c['old_panel']:c['missing_from_legacy_object'] for c in catalog if c['missing_from_legacy_object']},
    blacklisted={c['major']+'/'+c['old_panel']:c['blacklisted_genes'] for c in catalog if c['blacklisted_genes']},
    note='Minimal gene gate is only a diagnostic of the historical algorithm, not a biological validity test. Alternative-candidate sample replication is not recomputed here. V7 labels are unchanged.')
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
