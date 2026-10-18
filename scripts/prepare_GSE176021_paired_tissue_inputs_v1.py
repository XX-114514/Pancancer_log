"""Freeze exact equal-material case/tissue values and all paired sensitivities."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

MODES={'primary':('score','full'),'all_capture_mean':('all_captured_mean_sensitivity','full'),'fixed15_stress_complement':('score','fixed_stress_complement'),'own_population_background':('own_population_background_score','full'),'unresidualized_weighted_RNA':('unresidualized_weighted_RNA_score','full'),'RNA_quality_only_background':('RNA_quality_only_background_score','full')}


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);native=root/'GSE176021_author_lineage_native_RNA_20261019_v1';na=root/'GSE176021_author_lineage_native_RNA_audit_20261019_v1/private_terminal.json';assert sha(na)=='6950e1c15fda41b93962cfd763348898802d8875f1eef6719ec395b25392730d';scope=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';sa=root/'GSE176021_author_lineage_RNA_scope_audit_20261019_v1/private_terminal.json';assert sha(sa)=='fe5d399430fb1ced7637bb17200ef2327c794ae808d4653a77cbd1238504b982'
    for src in [native,scope]:
        t=json.loads((src/'private_terminal.json').read_text())
        for f,v in t['files'].items():assert sha(src/f)==v['sha256']
    out=root/'GSE176021_paired_tissue_inputs_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_case_differences_or_true_tests.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'native_acceptance_sha256':sha(na),'scope_acceptance_sha256':sha(sa),'primary':'All6 previously frozen questions; exact shared10case units, equal-material means within case/tissue on same-material paired-view eligible inputs. Original discovery-case SDs, no new learned scales.','sensitivity_modes':MODES,'author_confirmed_capture':'Full10-case sensitivity missing one pair stays incomplete; matched9 C/B point comparisons stored separately, never changes primary denominator, no new P.','coactivity':'Complete scope/case/tissue table permits point/LOO-only contextual review after full numeric acceptance; no repeated source P or MPR association.'},indent=2)+'\n')
    m=pd.read_csv(native/'private_all_author_scope_programme_measurements.tsv',sep='\t');pairs=pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv',sep='\t',keep_default_na=False);pairs=pairs[pairs.both_views_eligible];meta=pairs[['scope','sampleID','case_block','tissue','reviewed_patient']];m=m.merge(meta,on=['scope','sampleID'],validate='many_to_one');rows=[]
    for mode,(field,variant) in MODES.items():
        d=m[m.variant.eq(variant)].copy();d['value']=d[field];summary=d.groupby(['scope','case_block','tissue','program','view'],sort=True).agg(value=('value','mean'),materials=('value','size')).reset_index();summary['mode']=mode;rows.append(summary)
    case=pd.concat(rows,ignore_index=True);assert not case.duplicated(['scope','case_block','tissue','program','view','mode']).any();case.to_csv(out/'private_all_scope_case_tissue_values.tsv',sep='\t',index=False);ci=case.set_index(['scope','case_block','tissue','program','view','mode'])
    paired=pd.read_csv(scope/'private_all_paired_case_positions.tsv',sep='\t',keep_default_na=False);primary=paired[paired.scope.eq('all_author_lineage')].case_block.to_numpy(str);assert len(primary)==10
    old=root/'immune_RNA_MCP_new_source_family_inputs_20261019_v2/private_all6_original_discovery_case_scales.tsv';scales=pd.read_csv(old,sep='\t').set_index(['view','program']).discovery_case_sd
    registry=pd.read_csv(scope/'private_all6_paired_tissue_registered_tests.tsv',sep='\t');sensitivity=[];delta_records=[]
    def values(sc,cases,program,mode):
        R=[];E=[]
        for c in cases:
            for view,arr in [('R_FOXP3_capture',R),('E_effector_capture',E)]:
                index=[(sc,c,t,program,view,mode) for t in ['adjacent normal lung','tumor']]
                if all(k in ci.index for k in index):
                    n,t=[float(ci.loc[k,'value']) for k in index];arr.append((t-n)/scales.loc[(view,program)])
                else:arr.append(np.nan)
        return np.asarray(R),np.asarray(E)
    for r in registry.itertuples(index=False):
        R,E=values('all_author_lineage',primary,r.program,'primary');assert np.isfinite(R).all() and np.isfinite(E).all();x=R if r.outcome=='RE_tissue_cochange' else np.zeros(len(primary));y=E if r.outcome!='R_tissue_shift' else R
        np.savez_compressed(out/('private_'+r.model_id+'_inputs.npz'),cases=primary,x=x,y=y)
        for mode in MODES:
            a,b=values('all_author_lineage',primary,r.program,mode);xx=a if r.outcome=='RE_tissue_cochange' else np.zeros(len(primary));yy=b if r.outcome!='R_tissue_shift' else a
            for c,v,w in zip(primary,xx,yy):sensitivity.append({'model_id':r.model_id,'mode':mode,'case_block':c,'x':v,'y':w,'primary_N':10,'no_new_P':True})
        a,b=values('author_confirmed_capture',primary,r.program,'primary');xx=a if r.outcome=='RE_tissue_cochange' else np.zeros(len(primary));yy=b if r.outcome!='R_tissue_shift' else a
        for c,v,w in zip(primary,xx,yy):sensitivity.append({'model_id':r.model_id,'mode':'author_confirmed_capture_full_primary_cases','case_block':c,'x':v,'y':w,'primary_N':10,'no_new_P':True})
        shared=sorted(set(primary)&set(paired[paired.scope.eq('author_confirmed_capture')].case_block));assert len(shared)==9
        for sc,label in [('all_author_lineage','primary_on_author_confirmed_matched_subset'),('author_confirmed_capture','author_confirmed_matched_subset')]:
            a,b=values(sc,shared,r.program,'primary');xx=a if r.outcome=='RE_tissue_cochange' else np.zeros(len(shared));yy=b if r.outcome!='R_tissue_shift' else a
            for c,v,w in zip(shared,xx,yy):sensitivity.append({'model_id':r.model_id,'mode':label,'case_block':c,'x':v,'y':w,'primary_N':10,'no_new_P':True})
    registry.to_csv(out/'private_all6_registered_primary_tests.tsv',sep='\t',index=False);pd.DataFrame(sensitivity).to_csv(out/'private_all_prespecified_case_sensitivity_inputs.tsv',sep='\t',index=False);scales.rename('discovery_case_sd').to_csv(out/'private_original_discovery_case_scales.tsv',sep='\t')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    t={'status':'COMPLETE_ALL6_PAIRED_INPUTS_ALL_SCOPE_CASE_CONTEXTS_AND_SENSITIVITY_INPUTS_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'paired_primary_case_units':10,'models':6,'all_context_values':len(case),'sensitivity_case_rows':len(sensitivity),'new_P_or_clinical_association':False}
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps({'status':t['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
