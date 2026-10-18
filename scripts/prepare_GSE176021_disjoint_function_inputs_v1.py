"""Effect-before contract and every case input for all 400 functional questions."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


VIEWS = ['R_FOXP3_capture','E_effector_capture']
MODES = ['ucell1000','ucell2000','mean_logCP10K','fixed_stress_complement','own_population_background','unresidualized_weighted_RNA','RNA_quality_only_background','response_quality_partial_two']


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root,4*1024**2)
    assert sha(root/'GSE176021_disjoint_function_native_RNA_audit_20261019_v1/private_terminal.json') == 'b07e68106a615cb1f32f3e33f7acac18b5624c4f3cb82aa0180857d768bad21f'
    out = root/'GSE176021_disjoint_function_inputs_20261019_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan = {'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'native_acceptance_sha256':'b07e68106a615cb1f32f3e33f7acac18b5624c4f3cb82aa0180857d768bad21f','all_questions':'All50 complementary references x2 programmes x2 predictor lineages x2 response lineages=400 complete slots. Original10 paired cases, equal material within case/tissue, frozen programme scale. Within/cross role explicit, no manual pathway or direction selection.','primary':'Two-sided Spearman between frozen primary MCP tumour-minus-normal change and response-lineage complementary UCell1500 change. 99999 whole-case permutations per estimable slot, exact integer doubled-centred rank numerators; (k+1)/(B+1), two-sided integer exceedance. Frozen labelled seeds. All400 BH/BY include invalid slots asP1 internally, labelled missing effects/q remain missing.','uncertainty':'2000 shared whole-case-pair bootstrap; unadjusted marginal95percentile CI only if>=1900 finite; complete10leave-one effects. 99percent family simultaneous MCupper via .01/400 one-sided binomial bound then BH. Candidate flag requires BH<.05,MCupperBH<.05,finiteCIexcluding0 and allfiniteLOOsamesign. This flag is not receptor binding or clinical validation.','sensitivities':MODES,'quality_sensitivity':'Only a descriptive partial Spearman point with BOTH response-lineage changes of mean standardized logRNA and detected genes. Intercept and two rank quality controls mandatory, full rank3; no automatic dropping. Predictor/response residual zero marked unavailable. Not all8QC adjustment atN10, and no new sensitivity P.','limits':'Single NSCLC source, tissue contrasts include normal reference and possible shared upstream/treatment/inflammation. R/E programme synchronized predictors are not independent mechanisms; self-lineage vs cross-lineage functional associations reported together. Predictor union excluded from every functional target; gene independence does not eliminate quality or composition effects. Original references not adopted formal cell states. Frozen programme fitting uncertainty not included.','closure':'Old6tissue/420source/394directional families remain closed. No MPR endpoint, no clinical label adoption or expression download.'}
    (out/'plan_before_all400_functional_associations.json').write_text(json.dumps(plan,indent=2)+'\n')
    scope = root/'GSE176021_author_lineage_RNA_scope_20261019_v2'
    inp = root/'GSE176021_paired_tissue_inputs_20261019_v1'
    cases = np.load(inp/'private_AP01_inputs.npz')['cases']
    pairs = pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv',sep='\t')
    pairs = pairs[pairs.scope.eq('all_author_lineage') & pairs.case_block.isin(cases) & pairs.tissue.isin(['tumor','adjacent normal lung']) & pairs.both_views_eligible]
    assert len(pairs) == 60
    f = pd.read_csv(root/'GSE176021_disjoint_function_native_RNA_20261019_v1/private_all6000_functional_pool_values.tsv',sep='\t')
    f = f.merge(pairs[['sampleID','case_block','tissue']],on='sampleID',validate='many_to_one')
    refs = f.reference.drop_duplicates().tolist()
    assert len(refs) == 50 and len(f) == 6000
    value = pd.read_csv(inp/'private_all_scope_case_tissue_values.tsv',sep='\t')
    value = value[value.scope.eq('all_author_lineage')].set_index(['case_block','tissue','program','view','mode']).value
    scales = pd.read_csv(inp/'private_original_discovery_case_scales.tsv',sep='\t').set_index(['view','program']).discovery_case_sd
    native = pd.read_csv(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_all_author_scope_programme_measurements.tsv',sep='\t')
    native = native[native.scope.eq('all_author_lineage') & native.sampleID.isin(pairs.sampleID) & native.program.eq('MCP1') & native.variant.eq('full')]
    native = native.merge(pairs[['sampleID','case_block','tissue']],on='sampleID',validate='many_to_one')
    assert len(native) == 120
    quality = native.groupby(['case_block','tissue','view'])[['cellQ','geneQ']].mean()
    def pred(pr,view,mode):
        return np.array([(value.loc[(c,'tumor',pr,view,mode)]-value.loc[(c,'adjacent normal lung',pr,view,mode)])/scales.loc[(view,pr)] for c in cases])
    funcs = f.groupby(['reference','view','case_block','tissue'])[['ucell1500','ucell1000','ucell2000','mean_logCP10K']].mean()
    def target(ref,view,field):
        return np.array([funcs.loc[(ref,view,c,'tumor'),field]-funcs.loc[(ref,view,c,'adjacent normal lung'),field] for c in cases])
    rows=[];X=[];Y=[];SX=[];SY=[];Q=[]
    for ref in refs:
        for pr in ['MCP1','MCP2']:
            for predictor in VIEWS:
                for response in VIEWS:
                    x,y = pred(pr,predictor,'primary'),target(ref,response,'ucell1500')
                    sx=[];sy=[]
                    for mode in MODES:
                        sx.append(pred(pr,predictor,{'fixed_stress_complement':'fixed15_stress_complement'}.get(mode,mode)) if mode in MODES[3:7] else x)
                        sy.append(target(ref,response,mode) if mode in MODES[:3] else y)
                    q = np.array([(quality.loc[(c,'tumor',response)]-quality.loc[(c,'adjacent normal lung',response)]).to_numpy(float) for c in cases])
                    rows.append({'model_id':'CF'+str(len(rows)+1).zfill(3),'reference':ref,'program':pr,'predictor_view':predictor,'response_view':response,'role':'WITHIN_LINEAGE' if predictor == response else 'CROSS_LINEAGE','primary_case_units':10,'primary_materials':60,'metadata_preeligible':True,'family_denominator':400})
                    X.append(x);Y.append(y);SX.append(sx);SY.append(sy);Q.append(q)
    assert len(rows) == 400
    for array in [X,Y,SX,SY,Q]:assert np.isfinite(array).all()
    pd.DataFrame(rows).to_csv(out/'private_all400_registered_functional_questions.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all400_case_inputs.npz',cases=cases,model_ids=np.array([r['model_id'] for r in rows]),modes=np.array(MODES),x=np.array(X),y=np.array(Y),sensitivity_x=np.array(SX),sensitivity_y=np.array(SY),response_quality=np.array(Q))
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()}
    receipt={'status':'COMPLETE_ALL400_FUNCTIONAL_INPUTS_AND3200_POINT_SENSITIVITIES_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'primary_questions':400,'sensitivity_positions':3200,'real_functional_associations_executed':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__ == '__main__':main()
