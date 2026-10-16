"""Private complete functional maps; aggregate descriptions are not new inference."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def heatmaps(path,arrays,labels,refs,titles,footnote):
    fig=plt.figure(figsize=(21,16));axes=[];cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9')
    for ypos,values,title in zip([.62,.13],arrays,titles):
        ax=fig.add_axes([.18,ypos,.74,.29]);axes.append(ax);im=ax.imshow(np.ma.masked_invalid(values),vmin=-1,vmax=1,cmap=cmap,aspect='auto',interpolation='nearest')
        ax.set_yticks(np.arange(len(labels)),labels,fontsize=7.5);ax.set_xticks(np.arange(len(refs)),[x.removeprefix('HALLMARK_').replace('_',' ') for x in refs],rotation=70,ha='right',rotation_mode='anchor',fontsize=7)
        ax.tick_params(length=0);ax.set_title(title,loc='left',fontsize=12,pad=12);ax.set_xlim(-.5,len(refs)-.5);ax.set_ylim(len(labels)-.5,-.5)
    cb=fig.add_axes([.94,.30,.012,.48]);fig.colorbar(im,cax=cb,label='Descriptive partial rank correlation')
    fig.text(.18,.985,'NSCLC heldout-function RNA maps | complete fixed scope',ha='left',va='top',fontsize=15)
    fig.text(.18,.02,footnote,ha='left',va='bottom',fontsize=9)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
    visible=[]
    for ax in axes:
        visible.extend([ax.title,ax._left_title])
        for tick in ax.xaxis.get_major_ticks():
            if ax.get_xlim()[0]<=tick.get_loc()<=ax.get_xlim()[1]:visible.append(tick.label1)
        for tick in ax.yaxis.get_major_ticks():
            lo,hi=sorted(ax.get_ylim())
            if lo<=tick.get_loc()<=hi:visible.append(tick.label1)
    for t in visible:
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=width and b.y1<=height,(t.get_text(),'outside actual figure')
    fig.savefig(path,dpi=130,facecolor='white');plt.close(fig)
    return dict(visible_text_bounds='PASS_RENDERED_IN_VIEW_AXIS_LABELS_AND_PANEL_TITLES',rows=len(labels),columns=len(refs),panels=len(arrays))


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);paired=r/'NSCLC_matched_function_coordination_20261016_v1';proof=paired/'private_independent_matched_function_audit.json';deadline=time.monotonic()+14400
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_matched_function_coordination_v1/launch_receipt.json').read_text())
    while not proof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL28600') and sha(paired/'private_terminal.json')==json.loads(proof.read_text())['producer_terminal_sha256']
    qcroot=r/'NSCLC_RNA_composition_sensitivity_20261016_v1';qp=qcroot/'private_independent_sensitivity_audit.json';assert json.loads(qp.read_text())['status'].startswith('PASS_ALL212')
    out=r/'NSCLC_complete_functional_maps_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'private_display_plan.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),paired_audit_sha256=sha(proof),QC_audit_sha256=sha(qp),scope='All24 signatures x50 functions before/after QC and all572 pairs x50 functions, aggregated by original context/MCP. Medians not independent replication. No new P/BH or old-family changes.'),indent=2)+'\n')
    s=pd.read_csv(qcroot/'private_all1200_functional_QC_summary.tsv',sep='\t');t=pd.read_csv(paired/'private_all28600_matched_function_descriptions.tsv',sep='\t');refs=sorted(s.reference.unique());labels=sorted(s.signature.unique());assert len(refs)==50 and len(labels)==24
    before=s.pivot(index='signature',columns='reference',values='old_median_same_subtypes').reindex(index=labels,columns=refs);after=s.pivot(index='signature',columns='reference',values='QC_median_same_subtypes').reindex(index=labels,columns=refs)
    assert before.shape==after.shape==(24,50) and before.notna().all().all() and after.notna().all().all()
    bound1=heatmaps(out/'private_all24_functional_before_after_QC.png',[before.to_numpy(),after.to_numpy()],labels,refs,['A. Original 12 clinical and own composition/quality controls','B. Same original cases + author mitochondrial/ribosomal RNA percentages'],
                    'All 50 references exclude BOTH partners\' learned genes. Same cases shared across subtypes.\nRNA-pattern descriptions; no new P/q, protein-activity, formal state or independent-validation claim.')
    aggregates=[]
    for (context,program,reference),f in t.groupby(['context','program','reference'],sort=True):
        v=f[f.function_AB_rho_QC20.notna()];aggregates.append(dict(context=context,program=program,parent=f'C{context}_{program}',reference=reference,planned_subtype_pairs=len(f),describable_subtype_pairs=len(v),minimum_cases=int(v.eligible_paired_cases.min()) if len(v) else -1,maximum_cases=int(v.eligible_paired_cases.max()) if len(v) else -1,
              median_function_AB_QC20=v.function_AB_rho_QC20.median(),median_function_AB_original16=v.function_AB_rho_original16.median(),median_function_AB_plus_response=v.function_AB_rho_QC20_plus_response.median(),median_function_AB_plus_both_programmes=v.function_AB_rho_QC20_plus_both_programmes.median(),median_Aprogramme_Afunction=v.Aprogramme_Afunction_rho.median(),median_Bprogramme_Bfunction=v.Bprogramme_Bfunction_rho.median(),median_target_mean_AB=v.target_mean_AB_rho_QC20.median(),median_nontarget_mean_AB=v.nontarget_mean_AB_rho_QC20.median(),shared_cases_not_independent_replications=True))
    a=pd.DataFrame(aggregates);a.to_csv(out/'private_all_parent_function_descriptions.tsv',sep='\t',index=False,na_rep='NA');parents=sorted(a.parent.unique());p1=a.pivot(index='parent',columns='reference',values='median_function_AB_QC20').reindex(index=parents,columns=refs);p2=a.pivot(index='parent',columns='reference',values='median_function_AB_plus_both_programmes').reindex(index=parents,columns=refs)
    bound2=heatmaps(out/'private_all_parent_matched_function_coordination.png',[p1.to_numpy(),p2.to_numpy()],parents,refs,['A. Same-case function A-B: all 20 clinical and both-end composition/quality controls','B. Same function A-B + both original programme ranks (descriptive conditioning)'],
                    'Complete fixed pair/function scope; subtype-pair medians share cases and are not independent replications.\nConditioning on programmes is not causal mediation. No new P/q, directed communication or baseline prediction.')
    np.savez_compressed(out/'private_display_arrays.npz',signature_labels=np.asarray(labels,str),references=np.asarray(refs,str),single_before=before.to_numpy(),single_after=after.to_numpy(),parent_labels=np.asarray(parents,str),paired_main=p1.to_numpy(),paired_conditioned=p2.to_numpy())
    # Full closed-key display reconstruction, separate from upstream scientific OLS audits.
    disk=pd.read_csv(out/'private_all_parent_function_descriptions.tsv',sep='\t');assert len(disk)==t[['context','program','reference']].drop_duplicates().shape[0]
    for row in disk.itertuples():
        f=t[t.context.eq(row.context)&t.program.eq(row.program)&t.reference.eq(row.reference)];v=f[f.function_AB_rho_QC20.notna()];assert row.planned_subtype_pairs==len(f) and row.describable_subtype_pairs==len(v)
        assert np.allclose([row.median_function_AB_QC20,row.median_function_AB_plus_both_programmes,row.median_Aprogramme_Afunction,row.median_Bprogramme_Bfunction],[v.function_AB_rho_QC20.median(),v.function_AB_rho_QC20_plus_both_programmes.median(),v.Aprogramme_Afunction_rho.median(),v.Bprogramme_Bfunction_rho.median()],equal_nan=True)
    report=['# Complete NSCLC functional evidence maps (private)','',
            'All single-lineage and matched-pair functional slots are retained. All effects are descriptive, with no new P/BH. These shared-case medians are not meta-analysis or independent replication.',
            'Gene-excluded RNA patterns are not direct communication, measured protein activity, formal state labels, or baseline prediction. Author case identity is not fully reviewed.','',
            '## Complete maps','- private_all24_functional_before_after_QC.png','- private_all_parent_matched_function_coordination.png','',
            '## Full tables','- Existing all10600 QC sensitivities and all212 programme-composition descriptions remain in their immutable version directory.','- Existing all28600 matched-function rows remain in their immutable version directory.','- private_all_parent_function_descriptions.tsv retains every parent/function summary.','',
            '## Subsequent evidence','Connect only identical context/MCP/author endpoints to old clinical and coordination records; a single-end or generic functional pattern cannot replace failed paired clinical increment or absent direct LR evidence. Prioritize appropriate independent clinical/spatial validation after explicit measurement and case gates.']
    (out/'PRIVATE_COMPLETE_FUNCTIONAL_REVIEW.md').write_text('\n'.join(report)+'\n')
    terminal=dict(status='COMPLETE_PRIVATE_ALL24_AND_ALL_PARENT_FUNCTION_MAPS_ALL_DATA_CONSISTENCY_AND_RENDERED_LABEL_BOUNDS_PASS_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),single_map=bound1,paired_map=bound2,files={p.name:sha(p) for p in out.iterdir() if p.is_file()},no_new_P_BH=True,actual_visual_review=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
