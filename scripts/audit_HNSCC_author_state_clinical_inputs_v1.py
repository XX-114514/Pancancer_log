"""Rebuild clinical designs from original author counts using pandas rank and QR."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'HNSCC_author_state_clinical_preflight_20261016_v1';terminal=json.loads((pre/'private_terminal.json').read_text());plan=json.loads((pre/'plan_before_author_state_clinical_effects.json').read_text())
    for name,h in terminal['files'].items():assert sha(pre/name)==h
    source=r/'HNSCC_author_state_count_layer_preflight_20261016_v1';defs=pd.read_csv(source/'private_all35_author_immune_state_fraction_definitions.tsv',sep='\t').set_index('feature_id');raw=pd.read_csv(r/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t',dtype={'patient_registry':str,'library_id':str});clinical=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv';assert sha(clinical)==plan['clinical_sha256'];c=pd.read_csv(clinical,sep='\t',keep_default_na=False,dtype={'patient_id':str}).set_index('patient_id');scope=pd.read_csv(pre/'private_all276_state_clinical_scope.tsv',sep='\t',keep_default_na=False);model=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t').set_index('model_id');assert len(scope)==276 and scope.slot.is_unique and scope.family.value_counts().to_dict()==plan['scope'];measures=pd.read_csv(pre/'private_all_author_state_paired_measurements.tsv',sep='\t',dtype={'patient_id':str}).set_index(['feature_id','patient_id']);values={}
    for feature,row in defs.iterrows():
        labels=defs[defs.parent.eq(row.parent)].author_state.tolist();count={}
        for person,g in raw.groupby('patient_registry'):
            before=g[g.timepoint_registry.eq('pre')];after=g[g.timepoint_registry.eq('post')]
            if len(before)!=1 or len(after)!=1:continue
            a=before.iloc[0];b=after.iloc[0];na=sum(int(a[x]) for x in labels);nb=sum(int(b[x]) for x in labels)
            if min(na,nb)<20:continue
            baseline=int(a[row.author_state])/na;post=int(b[row.author_state])/nb;count[person]=(baseline,np.log1p(na),np.log1p(nb),post-baseline)
            m=measures.loc[(feature,person)];assert m.pre_library==a.library_id and m.post_library==b.library_id and m.pre_parent_cells==na and m.post_parent_cells==nb and m.pre_state_cells==int(a[row.author_state]) and m.post_state_cells==int(b[row.author_state]);assert np.allclose([m.baseline_fraction,m.post_fraction,m.delta_fraction],[baseline,post,post-baseline],rtol=1e-12,atol=1e-12)
        values[feature]=count
    assert len(measures)==sum(len(x) for x in values.values());seen=set();expected_links=[]
    for row in scope.itertuples():
        features=[row.feature_A] if not row.feature_B else [row.feature_A,row.feature_B];people=sorted(set.intersection(*(set(values[x]) for x in features)));assert len(people)==row.clinical_complete_patients
        for feature,label in zip(features,['A','B']):
            cc=c.loc[people];tx=cc.treatment.map({'D':0.,'D+T':1.}).to_numpy(float);hpv=cc.HPV.str.lower().map({'negative':0.,'positive':1.}).to_numpy(float);v=np.array([values[feature][x] for x in people]);d=np.column_stack([np.ones(len(people)),v[:,:3],tx,hpv,v[:,3]]);y=cc.regression.to_numpy(float)
            for j in [1,2,3,6]:
                ranked=pd.Series(d[:,j]).rank(method='average').to_numpy();d[:,j]=(ranked-ranked.mean())/ranked.std(ddof=0)
            yy=pd.Series(y).rank(method='average').to_numpy();yy=(yy-yy.mean())/yy.std(ddof=0);strata=(2*tx+hpv).astype(int);identifier=getattr(row,label+'_model')
            if not identifier:
                s=np.linalg.svd(d,compute_uv=False);assert len(y)<20 or len(y)-7<10 or not np.isfinite(d).all() or np.linalg.matrix_rank(d)<7 or s[0]/s[-1]>=1e8;continue
            p=pre/'private_designs'/(identifier+'.npz');assert sha(p)==model.loc[identifier,'design_sha256'];a=np.load(p,allow_pickle=False);assert a['patients'].astype(str).tolist()==people and np.array_equal(a['strata'],strata);assert np.allclose(a['D'],d,rtol=1e-12,atol=1e-12) and np.allclose(a['Y'],yy,rtol=1e-12,atol=1e-12)
            payload=json.dumps(people,separators=(',',':')).encode()+d.astype('<f8').tobytes()+yy.astype('<f8').tobytes()+strata.astype('<i8').tobytes();assert hashlib.sha256(payload).hexdigest()==identifier;assert model.loc[identifier,'n']==len(y) and model.loc[identifier,'residual_df']==len(y)-7
            q,rr=np.linalg.qr(d,mode='reduced');assert np.linalg.matrix_rank(rr)==7 and len(y)>=20 and len(y)-7>=10;assert np.isclose(model.loc[identifier,'max_leverage'],np.square(q).sum(axis=1).max(),rtol=1e-10,atol=1e-12);seen.add(identifier);expected_links.append(dict(slot=row.slot,feature_id=feature,model_id=identifier))
    links=pd.read_csv(pre/'private_all_slot_endpoint_model_links.tsv',sep='\t');pd.testing.assert_frame_equal(links.sort_values(['slot','feature_id']).reset_index(drop=True),pd.DataFrame(expected_links).sort_values(['slot','feature_id']).reset_index(drop=True),check_dtype=False);assert seen==set(model.index)
    proof=dict(status='PASS_ALL276_STATE_CLINICAL_SCOPE_ORIGINAL_COUNTS_EXACT_CASES_FULL7_COLUMN_DESIGNS_AND_ALL_MODELS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(pre/'private_terminal.json'),unique_models_audited=len(seen),logical_slots_audited=276,new_observed_clinical_effects=False);p=pre/'private_independent_state_clinical_input_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
