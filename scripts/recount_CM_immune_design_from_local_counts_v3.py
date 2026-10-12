"""Recount all v1/v2 aggregate designs from separately staged local counts."""
import argparse,csv,json,hashlib,itertools
from pathlib import Path
import numpy as np,pandas as pd

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workspace',required=True);a=ap.parse_args();w=Path(a.workspace);base=w/'annotation_execution_20261010';o=w/'Pancancer_log/inventories/analysis/2026-10-12/CM_immune_recovery_v1'
    cat={x['subCluster']:x['majorCluster'] for x in csv.DictReader((base/'downstream_review_v2/number_to_annotation_dictionary.tsv').open(),delimiter='\t')};names=sorted(cat)
    metap=base/'GSE243013_preflight_v1/immune_same_name_family_source_review_v1/family_clinical_context_review_v1/review_v2/all1290_family_context_scope.tsv'
    meta=list(csv.DictReader(metap.open(),delimiter='\t'));mi={x['sampleID']:x for x in meta};counts={s:{n:0 for n in names} for s in mi}
    cp=base/'downstream_review_v2/sample_major_subcluster_counts.tsv'
    for x in csv.DictReader(cp.open(),delimiter='\t'):
        if x['sampleID'] in counts:counts[x['sampleID']][x['subCluster']]+=int(x['n_cells'])
    assert all(sum(counts[s].values())==int(mi[s]['n_cells']) for s in mi)
    coretable=pd.read_csv(o/'fixed_CM_graph_core.tsv',sep='\t');cores={c:sorted(d.node) for c,d in coretable.groupby('cm')};weights={(x.cm,x.node):x.core_weight for x in coretable.itertuples()}
    immune={'B cell','T cell','T/NK lymphoid','Lymphoid cell','Dendritic cell','Macrophage','Myeloid','Monocyte','Neutrophil','Mast cell','Hematopoietic cell'}
    comp={'B cell':'B','T cell':'T','Dendritic cell':'myeloid','Macrophage':'myeloid','Myeloid':'myeloid','Monocyte':'myeloid','Neutrophil':'myeloid','Mast cell':'myeloid','Endothelial cell':'vascular','Fibroblast':'stromal','Mural cell':'stromal','Mesenchymal cell':'stromal','Epithelial cell':'epithelial','Melanocytic cell':'melanocytic','Glial cell':'glial'}
    def parent(m):
        if m in {'Dendritic cell','Macrophage','Myeloid','Monocyte','Neutrophil','Mast cell'}:return 'myeloid_analysis'
        if m in {'T cell','T/NK lymphoid','Lymphoid cell'}:return 'lymphoid_analysis'
        return m
    v1=pd.read_csv(o/'all_CM_cohort_results.tsv',sep='\t');v2=pd.read_csv(o/'rewiring_v2/all432_rewiring_sensitivities.tsv',sep='\t');framecache={};annotated=[];maxdiff=0;checks=0
    for key,d in v1.groupby('cohort'):
        gg,ca,cl,ph=key.split('|');samples=[x['sampleID'] for x in meta if x['gse_id']==gg and x['cancerType']==ca and x['sample_class']==cl and x['treatment_phase']==ph]
        assert len(samples)==d.libraries.iloc[0]
        allowed={n for n in names if gg!='GSE179994' or cat[n] in {'T cell','T/NK lymphoid','Lymphoid cell'}}
        for x in d.itertuples():
            core=cores[x.cm];rows=[]
            for s in samples:
                present=[n for n in core if n in allowed and counts[s][n]>=20];cw=sum(weights[x.cm,n] for n in present);maj=len({cat[n] for n in present});co=len({comp[cat[n]] for n in present if cat[n] in comp});minmaj=min(3,len({cat[n] for n in core}));minco=min(2,len({comp[cat[n]] for n in core if cat[n] in comp}))
                supported=len(present)>=int(np.ceil(len(core)/2)) and maj>=minmaj and co>=minco and cw>=.6;rows.append((len(present),cw,maj,supported,mi[s]['analysis_block']))
            expected=[np.median([z[0] for z in rows]),np.median([z[1] for z in rows]),np.median([z[2] for z in rows]),sum(z[3] for z in rows),len({z[4] for z in rows if z[3]})];actual=[x.median_core_ge20,x.median_core_weight,x.median_major,x.supported_libraries,x.supported_blocks];diff=max(abs(np.array(expected)-actual));maxdiff=max(maxdiff,float(diff));assert diff<1e-12;checks+=len(expected)
        for mode in ['legacy_major','parent_analysis']:
            group={n:parent(cat[n]) if mode=='parent_analysis' else cat[n] for n in names};groups=sorted(set(group.values()));tot={s:{g:sum(counts[s][n] for n in names if group[n]==g) for g in groups} for s in samples};op={g:np.mean([tot[s][g]>=20 for s in samples])>=.5 for g in groups};block=sorted({mi[s]['analysis_block'] for s in samples})
            values={};number={}
            for k in block:
                ss=[s for s in samples if mi[s]['analysis_block']==k];values[k]={};number[k]={}
                for n in names:
                    vals=[counts[s][n]/tot[s][group[n]] for s in ss if tot[s][group[n]]>0]
                    values[k][n]=np.mean(vals) if vals and n in allowed else np.nan;number[k][n]=np.mean([counts[s][n] for s in ss])
            framecache[(key,mode)]=(allowed,op,group,block,values,number)
    for row in v2.to_dict('records'):
        core=[n for n in cores[row['cm']] if row['view']=='whole_core' or cat[n] in immune];allowed,op,group,block,values,number=framecache[(row['cohort'],row['denominator'])];measured=[n for n in core if n in allowed and op[group[n]]];complete=[k for k in block if measured and all(np.isfinite(values[k][n]) for n in measured)]
        varying=[n for n in measured if len(complete)>=10 and np.mean([number[k][n]>=row['node_cell_cut'] for k in complete])>=.5 and sum(number[k][n]>=row['node_cell_cut'] for k in complete)>=10 and np.ptp([values[k][n] for k in complete])>1e-12]
        pairs=[(s,t) for s,t in itertools.combinations(varying,2) if group[s]!=group[t]];constant=sum(np.ptp([values[k][n] for k in complete])<=1e-12 for n in measured) if complete else 0
        assert [len(complete),len(varying),len(core)-len(measured),len(pairs),constant]==[row['complete_blocks'],row['variable_members'],row['unmeasured_members'],row['candidate_pairs'],row['constant_members']];checks+=5
        row['variable_nodes']=';'.join(varying);row['reestimated_pairs']=';'.join(s+'~'+t for s,t in pairs)
        if row['status']=='estimated_posthoc_diagnostic':
            used=sorted({n for pair in pairs for n in pair});V=np.array([[values[k][n] for n in used] for k in complete]);pi=[(used.index(s),used.index(t)) for s,t in pairs]
            ident=hashlib.sha256((row['cohort']+'|'+row['denominator']+'|'+';'.join(used)+'|'+str(pi)).encode()+V.tobytes()).hexdigest();assert ident==row['statistical_key']
            Z=pd.DataFrame(V).rank(method='average').to_numpy();Z-=Z.mean(axis=0);Z/=np.sqrt((Z*Z).sum(axis=0));rho=np.mean([Z[:,i]@Z[:,j] for i,j in pi]);diff=abs(rho-row['mean_rho']);maxdiff=max(maxdiff,float(diff));assert diff<1e-12;checks+=2
        annotated.append(row)
    pd.DataFrame(annotated).fillna('NA').replace('', 'NA').to_csv(o/'rewiring_v2/all432_rewiring_sensitivities_annotated_v3.tsv',sep='\t',index=False)
    summary=dict(status='PASS_SEPARATELY_STAGED_COUNTS_ALL54_CAPTURE_ALL432_DESIGNS_ALL120_OBSERVED_STATISTICS',numeric_checks=checks,max_absolute_difference=maxdiff,source_counts_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),clinical_context_sha256=hashlib.sha256(metap.read_bytes()).hexdigest(),independent_input_parsing=True,all_null_draws_replayed=False,raw_cell_matrix_replayed=False,independent_biological_validation=False)
    (o/'independent_design_recount_v3.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
