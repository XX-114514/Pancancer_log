"""Post-result signed endpoint-expression partition; no new P values."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_frozen92_axis_io_v2 import prepare


def score(L,R):
    return np.log1p(np.asarray(L,np.float32).astype(float)*np.asarray(R,np.float32).astype(float)).astype(np.float32).astype(float)


def parts(L,R,l0,r0):
    baseline=score(l0,r0);joint=score(L,R);left=score(L,r0);right=score(l0,R)
    return .5*((left-baseline)+(joint-right)),.5*((right-baseline)+(joint-left)),baseline


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_all92_endpoint_expression_decomposition_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);inp=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';main=root/'NSCLC_frozen92_axis_source_consistency_20261019_v1';proof=root/'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json';assert sha(proof)=='79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7'
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    save('plan_before_complete_endpoint_decomposition.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='Post-result descriptive all92 axes, not only nominal/supported candidate. Preserve all original P/BH/BY and failures, no new inferential tests or endpoints.',reference='Peraxis ligand and receptor component-minimum medians over all129 endpoint-measurable original discovery records, no response conditioning. f(L,R)=original float32 log1p product. Two-order Shapley terms average each marginal expression change at opposite endpoint reference and observed value. Terms sum exactly f(L,R)-constant.',clinical='Divide parts by original frozen product-score SD, project each through original full clinical design. Intercept cancels reference constant; signed clinical parts must sum original adjusted score effect. Direct least-squares target coefficient separately verifies each part. Signed shares are expression attribution, not causal contribution, receptor activation or proof ligand effect.',classification='Absolute receiver share>=.8 / sender>=.8 are post-result descriptive flags only. Shared signs/cancellation reported. All output private, no names/weights/individuals exported.'))
    rng=np.random.default_rng(2026101940);L=np.r_[0.,rng.lognormal(size=49)];R=np.r_[0.,rng.lognormal(size=49)];left,right,b=parts(L,R,np.median(L),np.median(R));np.testing.assert_allclose(left+right,score(L,R)-b,atol=1e-12,rtol=1e-12);save('private_synthetic_zero_and_complete_partition_gate.json',dict(status='PASS_EXACT_TWO_ENDPOINT_PARTITION_WITH_TRUE_ZERO'))
    a=dict(np.load(inp/'private_all129_all92_native_axis_and_design_inputs.npz'));vocab=json.loads((root/'NSCLC_communication_tensor_family_20261019_v1/private_selected_LR_vocabulary.json').read_text());raw=root/'NSCLC_full_native_gene_units_20261016_v1';means=np.load(raw/'mean_logCP10K.npy',mmap_mode='r');result=pd.read_csv(main/'private_all92_source_consistency_results.tsv',sep='\t').set_index('axis_id');rows=[];maximum=0.
    for j,axis in enumerate(a['axis_ids']):
        d=vocab[int(a['LR_index'][j])];s=int(a['sender_index'][j]);r=int(a['receiver_index'][j]);take=a['endpoint_measurable'][:,j];lv=np.zeros(129);rv=np.zeros(129)
        for i in np.flatnonzero(take):
            us=int(a['original_unit_indices'][i,s]);ur=int(a['original_unit_indices'][i,r]);lv[i]=np.float32(min(float(means[us,col]) for col in d['native_ligand_columns']));rv[i]=np.float32(min(float(means[ur,col]) for col in d['native_receptor_columns']))
        np.testing.assert_allclose(score(lv[take],rv[take]),a['raw_axis_score'][take,j],atol=0,rtol=0);l0=float(np.median(lv[take]));r0=float(np.median(rv[take]));pL,pR,baseline=parts(lv,rv,l0,r0);np.testing.assert_allclose((pL+pR)[take],a['raw_axis_score'][take,j]-baseline,atol=1e-12,rtol=1e-12);m,g,status,valid=prepare(a,j);row=dict(axis_id=axis,status=status,source='GSE243013',main_effect=result.loc[axis].effect,ligand_expression_part=np.nan,receptor_expression_part=np.nan,reconstructed_main_effect=np.nan,ligand_signed_share=np.nan,receptor_signed_share=np.nan,interpretation='UNESTIMABLE')
        if g is not None and result.loc[axis].status=='ESTIMABLE':
            sd=a['blind_score_mean_SD'][j,1];vectors=[pL[valid]/sd,pR[valid]/sd];beta=np.array([g['row']@y for y in vectors]);direct=np.array([np.linalg.lstsq(m['D'],y,rcond=1e-12)[0][-1] for y in vectors]);np.testing.assert_allclose(beta,direct,atol=1e-9,rtol=1e-9);np.testing.assert_allclose(beta.sum(),row['main_effect'],atol=1e-9,rtol=1e-9);maximum=max(maximum,float(abs(beta.sum()-row['main_effect'])));mass=np.abs(beta)/np.abs(beta).sum();category='RECEIVER_EXPRESSION_DOMINATES_SIGNED_PARTITION' if mass[1]>=.8 else 'SENDER_EXPRESSION_DOMINATES_SIGNED_PARTITION' if mass[0]>=.8 else 'BOTH_ENDPOINT_EXPRESSIONS_CONTRIBUTE_SIGNED';row.update(ligand_expression_part=float(beta[0]),receptor_expression_part=float(beta[1]),reconstructed_main_effect=float(beta.sum()),ligand_signed_share=float(beta[0]/row['main_effect']) if row['main_effect'] else np.nan,receptor_signed_share=float(beta[1]/row['main_effect']) if row['main_effect'] else np.nan,interpretation=category)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out/'private_all92_endpoint_expression_clinical_contributions.tsv',sep='\t',index=False);files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(d['bytes'] for d in files.values())<2*1024**2;receipt=dict(status='COMPLETE_ALL92_POST_RESULT_ENDPOINT_EXPRESSION_PARTITIONS_NO_NEW_P_VALUES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),original_result_audit_sha256=sha(proof),files=files,maximum_main_effect_reconstruction_error=maximum,original_P_BH_BY_changed=False,new_P_values=False,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
