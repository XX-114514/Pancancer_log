"""Literature-guided, multi-axis evidence reports. Never writes cell labels.

No argmax across identities/states, no contamination-string veto, no synthetic
probabilities. Expression, differential evidence, measurement and context remain
separate. The Lodi-like DEG profile is a reporting filter, not an identity gate.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def validate_registry(registry):
    panels=registry['panels']; byid={p['id']:p for p in panels}
    if len(byid)!=len(panels):raise ValueError('Duplicate panel IDs')
    refs={r['id'] for r in registry['references']}
    for p in panels:
        if p['auto_final_label'] or 'contamination' in p['id'].lower():raise ValueError('Unsafe panel semantics')
        if not set(p['parent_any'])<=set(byid)|{'*'}:raise ValueError('Missing parent')
        if not set(p['reference_ids'])<=refs:raise ValueError('Missing reference')
        if set(p['positive_genes'])&set(p['expected_lower_genes']):raise ValueError('Contradictory signs')
        if not set(sum(p['core_groups'],[])+p['naming_genes'])<=set(p['positive_genes']):raise ValueError('Invalid core/naming genes')
    def visit(k,stack):
        if k in stack:raise ValueError('Parent cycle')
        for par in byid[k]['parent_any']:
            if par!='*':visit(par,stack+[k])
    for k in byid:visit(k,[])
    return byid

def resolve_symbols(symbols,aliases):
    resolved=[aliases.get(str(g),str(g)) for g in symbols]
    if len(resolved)!=len(set(resolved)):
        raise ValueError('Alias/stable-ID collision; adjudicate source features, do not sum or drop automatically')
    return resolved

def standardize_expression(df):
    required={'cluster','gene','n_cells','n_measured','pct_expressing_measured'}
    if not required<=set(df):raise ValueError('Missing expression columns: '+str(required-set(df)))
    out=df.copy();out['cluster']=out.cluster.astype(str)
    for col in ['n_cells','n_measured']:
        out[col]=pd.to_numeric(out[col],errors='raise')
        if (~np.isfinite(out[col])|(out[col]<0)|(out[col]!=np.floor(out[col]))).any():raise ValueError('Invalid count denominator')
    if (out.n_cells==0).any():raise ValueError('Empty clusters are not expression observations')
    if out.duplicated(['cluster','gene']).any():raise ValueError('Duplicate cluster-gene rows')
    if ((out.n_measured<0)|(out.n_measured>out.n_cells)).any():raise ValueError('Invalid measured denominator')
    if out.groupby('cluster').n_cells.nunique().max()>1:raise ValueError('Inconsistent cluster counts')
    pct=pd.to_numeric(out.pct_expressing_measured,errors='raise')
    if ((out.n_measured>0)&(~np.isfinite(pct)|(pct<0)|(pct>1))).any():raise ValueError('Invalid measured detection fraction')
    out.loc[out.n_measured==0,'pct_expressing_measured']=np.nan
    return out

def evaluate(registry,expression,de=None,contexts=None,de_complete=False):
    byid=validate_registry(registry);expression=standardize_expression(expression)
    # Alias resolution is performed within cluster; source counts are never merged.
    for _,g in expression.groupby('cluster',sort=False):resolve_symbols(g.gene,registry['aliases'])
    expression['gene']=[registry['aliases'].get(g,g) for g in expression.gene]
    de_map={}
    if de is not None:
        req={'cluster','gene','avg_log2FC','pct.1','pct.2','p_val_adj'}
        if not req<=set(de):raise ValueError('DEG must be Seurat RNA avg_log2FC, not legacy mean-log difference')
        de=de.copy();de['cluster']=de.cluster.astype(str)
        for _,g in de.groupby('cluster',sort=False):resolve_symbols(g.gene,registry['aliases'])
        de['gene']=[registry['aliases'].get(g,g) for g in de.gene]
        if de.duplicated(['cluster','gene']).any():raise ValueError('DEG contrasts must be separated before evaluation')
        # Project positive-reporting filter on Lodi-style DE output. The
        # additional adjusted-p gate is NOT a claim of exact author logic.
        # Seurat p_val_adj is retained as supplied, not renamed to BH/FDR.
        ok=(de['avg_log2FC']>=.25)&(de['pct.1']>=.10)&((de['pct.1']-de['pct.2'])>=.25)&(de['p_val_adj']<.01)
        de_map={str(k):set(g.gene) for k,g in de[ok].groupby('cluster')}
    contexts=contexts or {};rows=[]
    for cluster,g in expression.groupby('cluster',sort=False):
        ix=g.set_index('gene');present={x for x in ix.index if ix.at[x,'n_measured']>0}
        deg=de_map.get(cluster,set());ctx=contexts.get(cluster,{})
        for p in registry['panels']:
            genes=p['positive_genes'];observed=[x for x in genes if x in present]
            pct={x:float(ix.at[x,'pct_expressing_measured']) for x in observed}
            components=[]
            for component in p['core_groups']:
                measured=[x for x in component if x in present]
                components.append(dict(requested=component,measured=measured,
                   max_pct=max((pct[x] for x in measured),default=None),
                   positive_DEG=sorted(set(component)&deg)))
            complete_components=all(x['measured'] for x in components)
            naming_measured=[x for x in p['naming_genes'] if x in present]
            naming_deg=sorted(set(p['naming_genes'])&deg)
            missing=[x for x in genes if x not in present]
            partial=[x for x in observed if ix.at[x,'n_measured']<ix.at[x,'n_cells']]
            # A known old parent is only context. A reviewed identity is a
            # separate external field, never inferred from these scores.
            parents=p['parent_any']
            accepted=set(ctx.get('reviewed_parent_ids',[]))
            parent_state='not_required' if not parents or '*' in parents else ('reviewed_context_present' if accepted&set(parents) else 'parent_not_independently_reviewed')
            context_state='not_required' if not p['required_context'] else ('context_present' if p['required_context'] in ctx.get('organ_contexts',[]) else 'organ_context_missing')
            if not observed:queue='measurement_unavailable'
            elif not complete_components:queue='core_component_unmeasured'
            elif p['evidence_level']=='curated_reference_gap':queue='reference_scope_review'
            elif context_state=='organ_context_missing':queue='organ_context_review'
            elif parent_state=='parent_not_independently_reviewed':queue='parent_review_then_panel_review'
            else:queue='multi_evidence_review'
            if not p['naming_genes']:name_state='no_marker_in_name_required'
            elif not naming_measured:name_state='naming_marker_unmeasured'
            elif naming_deg:name_state='named_marker_has_declared_DEG_evidence'
            elif any(pct[x]>0 for x in naming_measured):name_state='detected_not_proven_enriched'
            else:name_state='measured_without_positive_RNA_evidence'
            lower=p['expected_lower_genes']
            rows.append(dict(cluster=cluster,panel=p['id'],axis=p['axis'],n_cells=int(g.n_cells.iloc[0]),
                measured_genes=','.join(observed),unmeasured_or_not_in_input=','.join(missing),partially_measured_genes=','.join(partial),
                core_components_json=json.dumps(components,ensure_ascii=False),
                min_component_max_detection=min((x['max_pct'] for x in components if x['max_pct'] is not None),default=np.nan) if complete_components else np.nan,
                component_detection_is_joint_cell_fraction=False,
                each_core_component_has_measured_gene=complete_components,
                complete_core_measurement=all(x in present and ix.at[x,'n_measured']==ix.at[x,'n_cells'] for x in sum(p['core_groups'],[])),
                positive_DEG_genes=','.join(sorted(set(genes)&deg)),
                DEG_input_status='complete_declared_contrast' if de is not None and de_complete else ('truncated_DEG_input' if de is not None else 'DEG_not_supplied'),
                parent_review=parent_state,organ_context=context_state,naming_status=name_state,
                expected_lower_review='requires_measured_matched_comparator' if lower else 'not_required',
                expected_lower_measured=','.join(x for x in lower if x in present),
                review_queue=queue,final_label_assigned=False))
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--registry',type=Path,required=True);ap.add_argument('--expression',type=Path,required=True)
    ap.add_argument('--de',type=Path);ap.add_argument('--contexts',type=Path);ap.add_argument('--de-complete',action='store_true')
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    r=json.loads(a.registry.read_text());expr=pd.read_csv(a.expression,sep='\t',dtype={'cluster':str,'gene':str})
    de=pd.read_csv(a.de,sep='\t',dtype={'cluster':str,'gene':str}) if a.de else None
    ctx=json.loads(a.contexts.read_text()) if a.contexts else None
    result=evaluate(r,expr,de,ctx,a.de_complete)
    result.replace('',pd.NA).to_csv(a.out/'all_panel_evidence.tsv',sep='\t',index=False,na_rep='NA')
    manifest={'registry_sha256':hashlib.sha256(a.registry.read_bytes()).hexdigest(),'expression_sha256':hashlib.sha256(a.expression.read_bytes()).hexdigest(),
      'n_clusters':int(expr.cluster.nunique()),'n_panels':len(r['panels']),'rows':len(result),'final_cell_label_changes':0,
      'method':'all-candidate multi-axis evidence; no winner or margin; manual literature-guided interpretation',
      'de_filter':'Seurat supplied p_val_adj<.01; avg_log2FC>=.25; pct.1>=.1; pct.1-pct.2>=.25; not an identity classifier',
      'de_input_complete':a.de_complete,'limitations':['Aggregate detection does not establish same-cell coexpression','No common UCell threshold or automatic negative-marker exclusion','Context/DEG gaps remain explicit']}
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))
if __name__=='__main__':main()
