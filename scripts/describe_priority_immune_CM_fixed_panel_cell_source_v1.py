"""Cellular specificity of ALL existing fixed bulk panels across ALL native T nodes.

Post-result descriptive bridge; no clinical refit, new P or formal state labels.
"""
import json
import os
from datetime import datetime,timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_native_MFA_interpretation_inputs_v1 import sha


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,4*1024**2)
    clinical=Path(os.environ['PANCANCER_PROJECT_ROOT'])/'Project_v5/runs/20261018_priority_immune_CM_bulk_v1/analysis_v2'
    source=clinical/'native_T_panel_expression_private.tsv';contract=clinical/'frozen_plan_before_effects.json'
    fixed=json.loads(contract.read_text());table=pd.read_csv(source,sep='\t');nodes=sorted(table.node.unique())
    assert len(nodes)==16 and set(table.panel)==set(fixed['panels']) and not table.duplicated(['node','panel','gene']).any()
    out=root/'priority_immune_CM_fixed_panel_cell_source_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        existing_native_panel_file_sha256=sha(source),original_panel_contract_sha256=sha(contract),
        scope='Every original fixed panel and all16 original native T nodes. Explicitly post-result cellular interpretation, not new clinical inference or discovery of formal states.',
        expression_profile='For each gene, nonnegative node meanCP10K divided by its sum across all16 nodes; equal-gene mean within the unchanged panel. Separate abundance-weighted profile uses meanCP10K times measured cells, each gene normalized then averaged.',
        missing='Every required gene/node must have measured cells>0 and finite mean; all-zero genes or missing nodes make entire panel unavailable, not shorter. Record measured cell coverage heterogeneity.',
        comparison='All pairwise panel cosines and all nodes, top1/top3 concentration and entropy effective nodes. No P, threshold optimization or selection of spatial source based on effect.',
        caveat='Pooled node RNA and potentially unequal source capture are not same-cell coexpression, DEG, exact bulk patient percentile assay, original CM projection, local adjacency or secretion. Original annotation/model/clinical results unchanged.')
    (out/'plan_before_all_panel_cell_source_descriptions.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[];summaries=[];profiles={}
    for panel,genes in fixed['panels'].items():
        d=table[table.panel.eq(panel)]
        assert set(d.gene)==set(genes) and len(d)==len(nodes)*len(genes)
        expression=d.pivot(index='node',columns='gene',values='mean_CP10K').reindex(index=nodes,columns=genes).to_numpy(float)
        measured=d.pivot(index='node',columns='gene',values='measured_cells').reindex(index=nodes,columns=genes).to_numpy(float)
        expected=table.groupby('node').measured_cells.max().reindex(nodes).to_numpy(float)
        fraction=measured/expected[:,None];assert np.isfinite(fraction).all() and (fraction>=0).all() and (fraction<=1).all()
        valid=np.isfinite(expression).all() and (expression>=0).all() and (measured>0).all()
        if valid:valid=bool((expression.sum(axis=0)>0).all())
        if not valid:
            summaries.append(dict(panel=panel,status='UNAVAILABLE_COMPLETE_FIXED_PANEL',required_genes=len(genes),nodes=len(nodes),minimum_relative_measured_cell_coverage=float(fraction.min())))
            continue
        a=(expression/expression.sum(axis=0)).mean(axis=1)
        weighted=expression*measured;b=(weighted/weighted.sum(axis=0)).mean(axis=1)
        assert abs(a.sum()-1)<1e-12 and abs(b.sum()-1)<1e-12 and (a>=0).all() and (b>=0).all()
        profiles[panel]=(a,b);order=np.argsort(-a)
        summaries.append(dict(panel=panel,status='DESCRIPTIVE_FIXED_NODE_RNA_PROFILE',required_genes=len(genes),nodes=len(nodes),
            top_node=nodes[order[0]],top1_expression_share=float(a[order[0]]),top3_expression_share=float(a[order[:3]].sum()),
            effective_expression_nodes=float(np.exp(-np.sum(a[a>0]*np.log(a[a>0])))),
            minimum_relative_measured_cell_coverage=float(fraction.min()),maximum_relative_measured_cell_coverage=float(fraction.max()),
            all_gene_node_coverage_equal_to_available_node_max=bool(np.array_equal(measured,expected[:,None]*np.ones_like(measured)))))
        for i,node in enumerate(nodes):rows.append(dict(panel=panel,node=node,expression_profile_share=float(a[i]),measured_cell_weighted_RNA_share=float(b[i])))
    pairs=[]
    for x,y in combinations(fixed['panels'],2):
        if x not in profiles or y not in profiles:pairs.append(dict(panelA=x,panelB=y,status='UNAVAILABLE_COMPLETE_PANEL'));continue
        item=dict(panelA=x,panelB=y,status='DESCRIPTIVE_NOT_SAME_CELL_OR_CLINICAL_SCORE')
        for k,name in [(0,'equal_node_expression_cosine'),(1,'measured_cell_weighted_cosine')]:
            a,b=profiles[x][k],profiles[y][k];item[name]=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
        pairs.append(item)
    pd.DataFrame(rows).to_csv(out/'private_all_native_node_fixed_panel_profiles.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'private_all_fixed_panel_cell_source_summary.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(out/'private_all_panel_pair_specificity_descriptions.tsv',sep='\t',index=False)
    files=[p for p in out.iterdir() if p.is_file()]
    t=dict(status='COMPLETE_ALL_FIXED_PANEL_CELL_SOURCE_DESCRIPTIONS_NO_NEW_P',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),panels=len(fixed['panels']),nodes=len(nodes),pair_positions=len(pairs),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},independent_raw_or_new_clinical_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],completed_UTC=t['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
