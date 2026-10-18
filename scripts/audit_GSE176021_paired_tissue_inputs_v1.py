"""Independently reconstruct every case/tissue mean, SD and paired input."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def close(a,b):assert np.allclose(a,b,rtol=0,atol=1e-10,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);src=root/'GSE176021_paired_tissue_inputs_20261019_v1';out=root/'GSE176021_paired_tissue_input_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);term=wait(root,out,'prepare_GSE176021_paired_tissue_inputs_v1',src/'private_terminal.json',{'COMPLETE_ALL6_PAIRED_INPUTS_ALL_SCOPE_CASE_CONTEXTS_AND_SENSITIVITY_INPUTS_AUDIT_PENDING'})
    for n,v in term['files'].items():assert sha(src/n)==v['sha256']
    assert sha(root/'GSE176021_author_lineage_native_RNA_audit_20261019_v1/private_terminal.json')=='6950e1c15fda41b93962cfd763348898802d8875f1eef6719ec395b25392730d'
    scope=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';pairs=pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv',sep='\t',keep_default_na=False);native=pd.read_csv(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_all_author_scope_programme_measurements.tsv',sep='\t');table=pd.read_csv(src/'private_all_scope_case_tissue_values.tsv',sep='\t');modes={'primary':('score','full'),'all_capture_mean':('all_captured_mean_sensitivity','full'),'fixed15_stress_complement':('score','fixed_stress_complement'),'own_population_background':('own_population_background_score','full'),'unresidualized_weighted_RNA':('unresidualized_weighted_RNA_score','full'),'RNA_quality_only_background':('RNA_quality_only_background_score','full')};lookup={}
    for sc,cc,tissue,program,view,mode in table[['scope','case_block','tissue','program','view','mode']].itertuples(index=False,name=None):
        ids=pairs[(pairs.scope==sc)&(pairs.case_block==cc)&(pairs.tissue==tissue)&pairs.both_views_eligible].sampleID;field,variant=modes[mode];use=native[native.scope.eq(sc)&native.sampleID.isin(ids)&native.program.eq(program)&native['view'].eq(view)&native.variant.eq(variant)];assert len(use)==len(ids) and len(use)>0;value=float(np.sum(use[field].to_numpy(float))/len(use));key=(sc,cc,tissue,program,view,mode);lookup[key]=value
    assert len(lookup)==len(table)
    for r in table.itertuples(index=False):key=(r.scope,r.case_block,r.tissue,r.program,r.view,r.mode);close(r.value,lookup[key]);assert r.materials==int(((pairs.scope==r.scope)&(pairs.case_block==r.case_block)&(pairs.tissue==r.tissue)&pairs.both_views_eligible).sum())
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')&official.partition.eq('discovery')];scales={key:float(np.sqrt(np.mean((g.score-g.score.mean())**2))) for key,g in official.groupby(['view','program'])};saved=pd.read_csv(src/'private_original_discovery_case_scales.tsv',sep='\t').set_index(['view','program']).discovery_case_sd
    for k,v in scales.items():close(v,saved.loc[k])
    reg=pd.read_csv(src/'private_all6_registered_primary_tests.tsv',sep='\t');expected=pd.read_csv(scope/'private_all6_paired_tissue_registered_tests.tsv',sep='\t');pd.testing.assert_frame_equal(reg,expected);paired=pd.read_csv(scope/'private_all_paired_case_positions.tsv',sep='\t');cases=paired[paired.scope.eq('all_author_lineage')].case_block.to_numpy(str);assert len(cases)==10
    def delta(sc,cc,pr,view,mode):
        keys=[(sc,cc,t,pr,view,mode) for t in ['adjacent normal lung','tumor']]
        return (lookup[keys[1]]-lookup[keys[0]])/scales[(view,pr)] if all(k in lookup for k in keys) else np.nan
    for r in reg.itertuples(index=False):
        a=np.load(src/('private_'+r.model_id+'_inputs.npz'));assert np.array_equal(a['cases'],cases);R=np.array([delta('all_author_lineage',c,r.program,'R_FOXP3_capture','primary') for c in cases]);E=np.array([delta('all_author_lineage',c,r.program,'E_effector_capture','primary') for c in cases]);close(a['x'],R if r.outcome=='RE_tissue_cochange' else np.zeros(10));close(a['y'],R if r.outcome=='R_tissue_shift' else E)
    sens=pd.read_csv(src/'private_all_prespecified_case_sensitivity_inputs.tsv',sep='\t');idx=reg.set_index('model_id');shared=sorted(set(cases)&set(paired[paired.scope.eq('author_confirmed_capture')].case_block))
    for r in sens.itertuples(index=False):
        m=idx.loc[r.model_id];mode=r.mode;sc='all_author_lineage';base=mode
        if mode in ['author_confirmed_capture_full_primary_cases','author_confirmed_matched_subset']:sc='author_confirmed_capture';base='primary'
        elif mode=='primary_on_author_confirmed_matched_subset':base='primary'
        xx=delta(sc,r.case_block,m.program,'R_FOXP3_capture',base) if m.outcome=='RE_tissue_cochange' else 0.;yy=delta(sc,r.case_block,m.program,'R_FOXP3_capture' if m.outcome=='R_tissue_shift' else 'E_effector_capture',base);close([r.x,r.y],[xx,yy]);assert r.no_new_P and r.primary_N==10
    for (_,mode),g in sens.groupby(['model_id','mode']):assert g.case_block.tolist()==(shared if mode in ['author_confirmed_matched_subset','primary_on_author_confirmed_matched_subset'] else cases.tolist())
    d={'status':'ACCEPTED_COMPLETE_CASE_TISSUE_VALUES_ORIGINAL_SCALES_ALL6_INPUTS_AND_ALL_SENSITIVITIES','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'all_context_values':len(table),'all_sensitivity_case_rows':len(sens),'primary_pairs':10,'matched_author_confirmed_pairs':len(shared),'no_new_P_or_clinical_association':True}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
