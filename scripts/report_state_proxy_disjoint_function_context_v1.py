"""Full-family private interpretation; only source aggregates/figures exported."""
import hashlib,json,os,re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path('/dev/shm/pancancer-state-disjoint-functions-20261020-v1/private_output');audit=Path('/dev/shm/pancancer-state-disjoint-functions-audit-20261020-v2/private_output/private_terminal.json')
    a=json.loads(audit.read_text());assert a['auditor_sha256']=='a753da5084b7432bee5fb8d5ce77778ecaf1c00af78169c2afa9a80b9a30760a';assert a['status']=='ACCEPTED_ALL203_NATIVE_SCORE_RECONSTRUCTIONS_ALL_IDENTITY_LS_AND10150_SOURCE1550_AGGREGATES' and a['producer_terminal_sha256']==sha(root/'private_terminal.json')
    t=json.loads((root/'private_terminal.json').read_text())
    for n in ['private_all10150_source_function_descriptions.tsv.gz','private_all1550_cross_source_function_descriptions.tsv','private_frozen_gene_disjoint_function_contract.json']:assert sha(root/n)==t['files'][n]['sha256']
    source=pd.read_csv(root/'private_all10150_source_function_descriptions.tsv.gz',sep='\t',float_precision='round_trip');combined=pd.read_csv(root/'private_all1550_cross_source_function_descriptions.tsv',sep='\t',float_precision='round_trip');c=json.loads((root/'private_frozen_gene_disjoint_function_contract.json').read_text())
    assert len(source)==10150 and len(combined)==1550 and source.gse_id.map(lambda s:bool(re.fullmatch(r'GSE\d+',s))).all()
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    rows=[];annotations=[]
    for gid,b in combined.groupby('proxy_group_id',sort=True):
        supported=b[b.descriptive_concordance_flag];strong=supported[supported.strong_reference_representation]
        desc=[x for x in c['parent_descriptors'] if x['proxy_group_id']==gid]
        r=dict(proxy_group_id=gid,major=b.major.iloc[0],all_function_positions=len(b),positions_ge3_estimable_sources=int(b.estimable_sources.ge(3).sum()),
               descriptively_concordant_functions=len(supported),strong_representation_concordant_functions=len(strong),
               original_foreign_RNA_flag=any(type(x['foreign_RNA_sensitive']) is bool and x['foreign_RNA_sensitive'] for x in desc),
               original_foreign_RNA_flag_unknown=any(type(x['foreign_RNA_sensitive']) is not bool for x in desc),formal_state_assigned=False)
        annotations.append(r)
    annotation=pd.DataFrame(annotations);assert len(annotation)==31;annotation.to_csv(out/'private_complete31_function_interpretation_and_dictionary_limits.tsv',sep='\t',index=False)
    for key,b in source.groupby(['gse_id','major'],sort=True):
        z=b[b.estimable_identity_strata.gt(0)];r=dict(zip(['public_GSE','major'],key));r.update(predefined_proxy_positions=int(b.proxy_group_id.nunique()),predefined_function_positions=len(b),estimable_function_positions=len(z),both_adjusted_rhos_ge_point2=int((z.adjusted_cell_rho.ge(.2)&z.adjusted_block_rho.ge(.2)).sum()),both_adjusted_rhos_le_minus_point2=int((z.adjusted_cell_rho.le(-.2)&z.adjusted_block_rho.le(-.2)).sum()),strong_reference_estimable_positions=int(z.strong_reference_representation.sum()),new_independent_patient_validation=False)
        for col in ['raw_rho','adjusted_cell_rho','adjusted_block_rho']:
            for label,q in [('q25',.25),('median',.5),('q75',.75)]:r[label+'_'+col]=float(z[col].quantile(q)) if len(z) else np.nan
        rows.append(r)
    refsummary=pd.DataFrame(rows);assert len(refsummary)==12;refsummary.to_csv(out/'source_native_function_coverage_and_effects.tsv',sep='\t',index=False)
    models={};folds={}
    questions={'BACKGROUND_ONLY':'pancancer-real-spatial-mv-background-20261020-v1','OTHER_RECEIVER_STATES_CONDITIONAL':'pancancer-real-spatial-mv-20261020-v1'}
    expected_terms={'BACKGROUND_ONLY':'6249069ad1f86f80bfefddb61d66dcc64d9d3a0f5c44d20226553ec8f71b9916','OTHER_RECEIVER_STATES_CONDITIONAL':'e187616ee043701294421c002444973803684be0a1fd34d4a050c4fa4d5e0bbe'}
    for question,slug in questions.items():
        p=Path('/dev/shm')/slug/'private_output';assert sha(p/'private_terminal.json')==expected_terms[question];mt=json.loads((p/'private_terminal.json').read_text());n='private_complete126_blocked_model_cases.tsv';assert sha(p/n)==mt['files'][n]['sha256'];models[question]=pd.read_csv(p/n,sep='\t',float_precision='round_trip')
        folds[question]={}
        for i in range(4):
            n='private_image_%02d_folds.tsv.gz'%i;assert sha(p/n)==mt['files'][n]['sha256'];folds[question][i]=pd.read_csv(p/n,sep='\t',float_precision='round_trip')
    physical=json.loads(Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output/private_frozen_matched_measurement_and_target_contract.json').read_text());ai=annotation.set_index('proxy_group_id');target=[]
    for question,model in models.items():
        for i,t in enumerate(physical['targets']):
            for m in t['proxy_mappings']:
                z=model[model.image_ordinal.eq(i)&model.target_proxy_group_id.eq(m['proxy_group_id'])];assert len(z)==2;okay=bool(z.status.eq('ESTIMABLE_BLOCKED_PREDICTION').all());f=folds[question][i];f=f[f.target_proxy_group_id.eq(m['proxy_group_id'])] if len(f) else f
                r=dict(public_GSE=t['dataset_id'],cancer=t['cancer'],major=m['major'],question=question,proxy_group_id=m['proxy_group_id'],image_ordinal=i,
                    paired_estimable=okay,both_readouts_gain_positive=bool(okay and z.gain_R2.gt(0).all()),all_six_folds_gain_positive=bool(okay and len(f)==6 and f.gain_R2.gt(0).all()),
                    function_descriptively_concordant=bool(ai.loc[m['proxy_group_id'],'descriptively_concordant_functions']>0),strong_function_descriptively_concordant=bool(ai.loc[m['proxy_group_id'],'strong_representation_concordant_functions']>0),
                    original_foreign_RNA_flag=bool(ai.loc[m['proxy_group_id'],'original_foreign_RNA_flag']),original_foreign_RNA_flag_unknown=bool(ai.loc[m['proxy_group_id'],'original_foreign_RNA_flag_unknown']))
                target.append(r)
    target=pd.DataFrame(target);assert len(target)==126
    target.to_csv(out/'private_all126_target_question_function_and_prediction_descriptions.tsv',sep='\t',index=False)
    summaries=[]
    for key,b in target.groupby(['public_GSE','cancer','major','question'],sort=True):
        r=dict(zip(['public_GSE','cancer','major','question'],key));r.update(target_positions=len(b),paired_estimable_positions=int(b.paired_estimable.sum()),function_descriptively_concordant_positions=int(b.function_descriptively_concordant.sum()),strong_function_descriptively_concordant_positions=int(b.strong_function_descriptively_concordant.sum()),
          function_and_both_readout_gain_positions=int((b.function_descriptively_concordant&b.both_readouts_gain_positive).sum()),strong_function_and_both_readout_gain_positions=int((b.strong_function_descriptively_concordant&b.both_readouts_gain_positive).sum()),
          function_and_all_six_fold_gain_positions=int((b.function_descriptively_concordant&b.all_six_folds_gain_positive).sum()),strong_function_and_all_six_fold_gain_positions=int((b.strong_function_descriptively_concordant&b.all_six_folds_gain_positive).sum()),reviewed_independent_patients=0)
        summaries.append(r)
    targetsummary=pd.DataFrame(summaries);assert len(targetsummary)==8;targetsummary.to_csv(out/'target_source_function_spatial_support.tsv',sep='\t',index=False)
    limits=target[target.question.eq('BACKGROUND_ONLY')].groupby(['public_GSE','cancer','major'],as_index=False).agg(target_positions=('proxy_group_id','size'),original_foreign_RNA_flag_positions=('original_foreign_RNA_flag','sum'),original_foreign_RNA_flag_unknown_positions=('original_foreign_RNA_flag_unknown','sum'))
    limits.to_csv(out/'target_source_original_dictionary_limits.tsv',sep='\t',index=False)
    identity=pd.read_csv(root/'private_all_identity_disjoint_function_descriptions.tsv.gz',sep='\t',float_precision='round_trip');reasons=identity.groupby(['gse_id','major','status'],sort=True).size().rename('identity_function_positions').reset_index().rename(columns={'gse_id':'public_GSE'});reasons.to_csv(out/'source_identity_function_reason_counts.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(16,9));pos=np.arange(len(refsummary));labels=refsummary.public_GSE+' | '+refsummary.major.replace({'Macrophage':'Mac','T cell':'T'})
    for col,offset,color,label in [('adjusted_cell_rho',-.1,'#347eae','Within-identity cell rho'),('adjusted_block_rho',.1,'#b46431','Residual block-mean rho')]:
        mid=refsummary['median_'+col].to_numpy();lo=refsummary['q25_'+col].to_numpy();hi=refsummary['q75_'+col].to_numpy();axes[0].errorbar(mid,pos+offset,xerr=np.vstack([mid-lo,hi-mid]),fmt='o',color=color,label=label,capsize=2)
    axes[0].set_yticks(pos,labels);axes[0].axvline(0,color='#777777',lw=.8);axes[0].set_xlim(-1,1);axes[0].set_xlabel('Across-function median and IQR (not CI)');axes[0].set_title('Reference sources | gene-disjoint function descriptions')
    order=[]
    for cancer,major in [('LUAD','T cell'),('LUAD','Macrophage'),('pleural_mesothelioma','T cell'),('pleural_mesothelioma','Macrophage')]:
        for question in questions:order.append((cancer,major,question))
    z=targetsummary.set_index(['cancer','major','question']).loc[order];yy=np.arange(8)
    axes[1].barh(yy,z.target_positions,color='#e2e8eb',label='Frozen target positions')
    axes[1].barh(yy,z.function_descriptively_concordant_positions,color='#80afbd',label='Any concordant function complement')
    axes[1].barh(yy,z.function_and_both_readout_gain_positions,color='#cf8b46',label='Function + positive gain in both readouts')
    axes[1].plot(z.function_and_all_six_fold_gain_positions,yy,'|',color='#222222',ms=10,label='Function + positive gain in all six folds')
    ylabels=[('LUAD' if x[0]=='LUAD' else 'Mesothelioma')+' | '+('T' if x[1]=='T cell' else 'Mac')+' | '+('Background' if x[2]=='BACKGROUND_ONLY' else '+ own states') for x in order]
    axes[1].set_yticks(yy,ylabels);axes[1].set_xlabel('Correlated measurement positions; not distinct states');axes[1].set_title('GSE299786 | function context and full spatial family')
    for ax in axes:ax.invert_yaxis()
    h0,l0=axes[0].get_legend_handles_labels();h1,l1=axes[1].get_legend_handles_labels();fig.legend(h0+h1,l0+l1,loc='upper center',bbox_to_anchor=(.54,.925),ncol=2,frameon=False,fontsize=9)
    fig.suptitle('All measured proxies: independent-gene functional context',fontsize=15,y=.98)
    fig.text(.025,.035,'Fixed 31 proxies x 50 references; complete source and failed positions retained. Descriptive concordance has no P/FDR or formal state/causal meaning.\nRNA quality and identity strata controlled; the reference cells were also used for proxy qualification. This is not untouched external validation.\nSpatial T labels are broad author lineages, not confirmed gamma-delta identities. Reference blocks inherit prior units; spatial inputs have no reviewed patient bridge.',fontsize=9)
    fig.subplots_adjust(left=.16,right=.98,top=.76,bottom=.17,wspace=.8)
    for ext in ['png','svg']:fig.savefig(out/('source_function_and_spatial_context.'+ext),dpi=180)
    plt.close(fig)
    terminal=dict(status='COMPLETE_ALL31_FUNCTION_CONTEXT_AND_SOURCE_ONLY_REPORT_PENDING_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),function_audit_sha256=sha(audit),formal_state_or_interaction_or_clinical_validation=False,scientific_outputs_public=False,
                  files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()})
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps({'status':terminal['status'],'file_count':len(terminal['files'])}))


if __name__=='__main__':main()
