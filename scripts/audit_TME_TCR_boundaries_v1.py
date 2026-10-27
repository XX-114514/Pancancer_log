"""Independent complete replay of the additional TCR-core proposal policy.

No new production policy or proposal builder is imported. Frozen storage and
vocabulary are inputs, not independent identity truth.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

FAM=['T_NK','B_Plasma','Myeloid','Mast','Endothelial','Fibroblast','Mural','Epithelial','Glial','Astrocyte','Neuron','Erythroid']
PANELS={
 'B_Plasma':'CD79A CD79B MS4A1 CD19 CD22 PAX5 MZB1 JCHAIN DERL3 TNFRSF17',
 'Myeloid':'CSF1R SPI1 LST1 LYZ AIF1 CLEC9A XCR1 CD1C FCER1A LILRA4 CLEC4C IL3RA TCF4 CSF3R FCGR3B CEACAM8 CXCR2 MMP8 CAMP',
 'Mast':'KIT TPSAB1 TPSB2 CPA3 HDC MS4A2',
 'Endothelial':'PECAM1 CDH5 CLDN5 VWF ESAM PROX1 FLT4 CCL21 PDPN',
 'Fibroblast':'COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA',
 'Mural':'RGS5 CSPG4 PDGFRB MCAM MYH11 CNN1',
 'Epithelial':'EPCAM CDH1 KRT7 KRT8 KRT18 KRT19 KRT5 KRT14 SFN',
 'Glial':'OLIG1 OLIG2 SOX10 MBP PLP1 MOG MOBP MAG',
 'Astrocyte':'GFAP AQP4 SLC1A3','Neuron':'SNAP25 SYT1 RBFOX3 STMN2','Erythroid':''}
REASONS=['not_baseline_unresolved','additional_T_RNA_core_reference_local_concordance','TCR_complex_or_chain_or_signaling_not_complete','competing_lineage_core_detected','old_competing_family_core_measurement_incomplete','reference_family_not_unique_T_NK','low_reference_not_T_or_unassigned','source_local_group_not_T_or_T_NK','without_fitted_source_group','zero_native_RNA']

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
    begin=time.monotonic();P=Path(load(root/'config/spec_v1.json')['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    cfg=load(P/'config/spec_v1.json');E=Path(cfg['native_evidence']);J=Path(cfg['joined']);F=Path(cfg['candidate_producer'])
    for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
    for stage in [E,Path(cfg['native_audit']),Path(cfg['frozen_policy_replay'])]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);assert sha(w['index']['path'])==w['index']['SHA256'];wi={r['source']:r for r in load(w['index']['path'])};base=Path(w['index']['path']).parent
    ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);ri={r['source']:r for r in load(ref/'outputs/full_archive_reference_index_v1.json')};models=load(ref/'outputs/reference_models_v1.json')
    sys.path.insert(0,str(F/'scripts'));from native_lineage_routing_v1 import reference_family
    lowmap=[]
    for label in models['immune_low']['labels']:
        s=label.lower();eligible=reference_family('immune_low',label)=='T_NK'
        isT=('nk cells' not in s and s!='transitional nk' and not s.startswith('ilc') and ('t cells' in s or 'thymocytes' in s or s in {'cd8a/a','cd8a/b(entry)','t(agonist)','treg(diff)'}))
        lowmap.append(eligible and isT)
    lowmap=np.asarray(lowmap,bool);lmroot=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);li={r['source']:r for r in load(lmroot/'outputs/full_local_identity_review_index_v1.json')}
    t=load(P/'provenance/terminal_v1.json');ip=P/'outputs/TCR_boundary_proposal_index_v1.json';assert sha(ip)==t['index_SHA256'];pi={r['source']:r for r in load(ip)}
    ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};assert set(pi)==set(wi)==set(ei)
    totals=collections.Counter();reasons=collections.Counter();controls=collections.Counter()
    for ds in sorted(wi):
        pe=pi[ds];mp=P/'outputs'/pe['metadata'];assert sha(mp)==pe['metadata_SHA256'];m=load(mp);pa=P/'outputs'/m['arrays'];assert sha(pa)==m['arrays_SHA256']
        with np.load(pa,allow_pickle=False) as f:out={k:f[k] for k in f.files}
        wr=wi[ds];wm=load(base/wr['metadata']);assert sha(base/wr['metadata'])==wr['metadata_SHA256'] and sha(base/wr['arrays'])==wr['arrays_SHA256']
        with np.load(base/wr['arrays'],allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_family_bitmask','reference_family_bitmask','consensus_major_code','annotation_group_code','parent_source_local_group_code','fit_pool_candidate_not_final_QC']}
        em=E/'outputs'/ei[ds]['metadata'];assert sha(em)==ei[ds]['metadata_SHA256'];e=load(em);c=e['native_context'];assert sha(c['path'])==c['SHA256']
        with np.load(c['path'],allow_pickle=False) as f:x=f['native_gene_values'];genes=e['context_genes'];mask=f['measured'];rna=f['native_total_RNA'];np.testing.assert_array_equal(q['native_row'],f['native_row'])
        for k in ['native_row','source_native_position']:np.testing.assert_array_equal(q[k],out[k])
        rr=ri[ds];rp=ref/'outputs'/rr['metadata'];ra=ref/'outputs'/rr['arrays'];assert sha(rp)==rr['metadata_sha256'] and sha(ra)==rr['arrays_sha256'];rm=load(rp)
        detected={g:(x[:,j]>0 if mask[j] else np.zeros(len(x),bool)) for j,g in enumerate(genes)};available={g:bool(mask[j]) for j,g in enumerate(genes)};restored=[]
        with np.load(ra,allow_pickle=False) as f:
            pos=q['source_native_position'];np.testing.assert_array_equal(q['native_row'],f['native_row'][pos]);low=f['immune_low__unique_code'][pos];lowT=np.zeros(len(low),bool);valid=low>=0;lowT[valid]=lowmap[low[valid]]
            required={g for s in PANELS.values() for g in s.split()}
            for g in sorted(required):
                if g in genes and not available[g] and g in rm['review_genes'] and rm['gene_measured'][rm['review_genes'].index(g)]:
                    k=rm['review_genes'].index(g);detected[g]=np.unpackbits(f['native_gene_detection_packbits'][pos,k//8:k//8+1],axis=1,bitorder='little')[:,k%8].astype(bool);available[g]=True;restored.append(g)
        assert restored==m['binary_only_alias_restoration']
        z=np.zeros(len(x),bool)
        def g(name):return detected.get(name,z)
        def number(names):return np.column_stack([g(v) for v in names.split()]).sum(axis=1)
        alpha=g('TRAC')|g('TRBC1')|g('TRBC2');gamma=g('TRDC')&(g('TRGC1')|g('TRGC2'))
        core=g('CD3D')&g('CD3E')&(alpha|gamma)&(number('LCK LAT ZAP70 BCL11B')>0)
        b=((g('CD79A')|g('CD79B'))&(number('MS4A1 CD19 CD22 PAX5')>0))|(g('CD79A')&g('CD79B'))|(g('MZB1')&(number('JCHAIN DERL3 TNFRSF17')>0))
        mye=((g('CSF1R')|g('SPI1'))&(number('CSF1R SPI1 LST1 LYZ AIF1')>2))|(g('CLEC9A')&g('XCR1'))|(g('CD1C')&g('FCER1A'))|((g('LILRA4')|g('CLEC4C'))&(number('LILRA4 CLEC4C IL3RA TCF4')>2))|((g('CSF3R')|g('FCGR3B'))&(number('CEACAM8 CXCR2 MMP8 CAMP')>0))
        ccore={'B_Plasma':b,'Myeloid':mye,'Mast':(number('TPSAB1 TPSB2 CPA3')>0)&(number(PANELS['Mast'])>1),
            'Endothelial':(number('PECAM1 CDH5 CLDN5 VWF ESAM')>2)|((g('PROX1')|g('FLT4'))&(number('PROX1 FLT4 CCL21 PDPN')>2)),
            'Fibroblast':number(PANELS['Fibroblast'])>2,'Mural':(g('RGS5')|g('CSPG4')|g('MYH11'))&(number(PANELS['Mural'])>2),
            'Epithelial':number(PANELS['Epithelial'])>2,'Glial':(number('OLIG1 OLIG2 SOX10')>0)&(number(PANELS['Glial'])>2),
            'Astrocyte':number(PANELS['Astrocyte'])>1,'Neuron':number(PANELS['Neuron'])>1,'Erythroid':(q['native_family_bitmask']&2048)>0}
        competitor=np.zeros(len(x),np.uint16);unknown=np.zeros(len(x),np.uint16);coverage={}
        for name,d in ccore.items():
            bit=1<<FAM.index(name);competitor[d]|=bit;complete=all(available.get(k,False) for k in PANELS[name].split());coverage[name]=complete
            if not complete:unknown[(q['native_family_bitmask']&bit)>0]|=bit
        assert coverage==m['competing_core_measurement_complete']
        lp=lmroot/'outputs'/li[ds]['output'];assert sha(lp)==li[ds]['output_sha256'];localset={r['source_local_group_code'] for r in load(lp)['groups'] if r['local_major_proposal'] in ['T','T_NK']};localT=np.isin(q['parent_source_local_group_code'],list(localset))
        assert m['inputs']=={'working_metadata_SHA256':sha(base/wr['metadata']),'working_arrays_SHA256':sha(base/wr['arrays']),'native_context_metadata_SHA256':sha(em),'reference_metadata_SHA256':sha(rp),'reference_arrays_SHA256':sha(ra),'local_review_SHA256':sha(lp)}
        major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];unresolved=major=='unresolved';code=np.zeros(len(x),np.uint8)
        # Reverse assignment gives the same explicit precedence, with a
        # different implementation from the production pending-mask loop.
        code[unresolved]=1
        for value,bad in [(7,~localT),(8,~q['fit_pool_candidate_not_final_QC']),(6,~lowT),(5,q['reference_family_bitmask']!=1),(4,unknown!=0),(3,competitor!=0),(2,~core),(9,rna<=0)]:code[unresolved&bad]=value
        proposed=code==1
        expected={'review_code':code,'proposed_T':proposed,'TCR_core':core,'alpha_chain_detected':alpha,'gamma_chain_detected':gamma,'competing_core_mask':competitor,'unassessed_old_competitor_mask':unknown,'reference_low_T':lowT,'source_local_T_context':localT,'baseline_unresolved':unresolved}
        for key,v in expected.items():np.testing.assert_array_equal(out[key],v)
        reason=collections.Counter(REASONS[int(v)] for v in code[unresolved]);assert dict(reason)==m['reason_counts_in_baseline_unresolved']
        affected={int(v['source_annotation_group_code']):v for v in m['affected_previous_exact_groups']};assert set(affected)==set(map(int,np.unique(q['annotation_group_code'][proposed])))
        for k,r in affected.items():
            old=q['annotation_group_code']==k;take=old&proposed;assert r['n_previous_exact_group']==int(old.sum()) and r['n_proposed_T']==int(take.sum()) and r['n_remainder']==int((old&~take).sum())
            for gene in ['CD3D','CD3E']:
                v=x[take,genes.index(gene)];s=r['own_proposed_subset_fixed_marker_expression_not_DEG'][gene];assert (v>0).all() and s['detection_fraction']==1
                np.testing.assert_allclose(s['mean_native_count'],v.sum()/len(v),rtol=1e-13)
                np.testing.assert_allclose(s['native_mean_CP10K'],sum(v*10000/rna[take])/len(v),rtol=1e-12)
        for label in ['T','NK','Myeloid','Epithelial','unresolved']:
            sel=major==label;v={'n_records':int(sel.sum()),'n_with_TCR_core':int((sel&core).sum()),'n_with_TCR_and_competing_core':int((sel&core&(competitor!=0)).sum())};assert v==m['baseline_label_control_counts_not_accuracy'][label]
            for key,value in v.items():controls[(label,key)]+=value
        assert m['n_scope_records']==len(x) and m['n_baseline_unresolved']==int(unresolved.sum()) and m['n_proposed_T']==int(proposed.sum())==pe['n_proposed_T']
        totals.update(n_sources=1,n_scope_records=len(x),n_baseline_unresolved=int(unresolved.sum()),n_proposed_T=int(proposed.sum()),n_affected_previous_groups=len(affected),n_split_previous_groups=sum(r['n_remainder']>0 for r in affected.values()));reasons.update(reason)
        save(root/'outputs'/(ds+'_independent_TCR_boundary_replay_v1.json'),{'source':ds,'n_scope_records':len(x),'all_policy_fields_group_subsets_and_fixed_marker_statistics_replayed':True,'n_proposed_T':int(proposed.sum())})
    assert dict(totals)==t['totals'];sp=P/'outputs/TCR_boundary_proposal_summary_v1.json';assert sha(sp)==t['summary_SHA256'];s=load(sp);assert s['totals']==dict(totals) and s['reason_counts_in_baseline_unresolved']==dict(reasons)
    for label,r in s['baseline_control_counts_not_accuracy'].items():
        for k,v in r.items():assert controls[(label,k)]==v
    report={'status':'COMPLETE_CURRENT_SCOPE_TCR_CORE_PROPOSALS_INDEPENDENTLY_REPLAYED_NOT_APPLIED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'totals':dict(totals),'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),'producer_index_SHA256':sha(ip),'producer_summary_SHA256':sha(sp),
        'all_policy_arrays_exactly_replayed':True,'new_production_policy_not_imported':True,'not_biological_accuracy_or_purity_certification':True,'current_consensus_QC_states_and_CM_unchanged':True,'wall_seconds':time.monotonic()-begin}
    op=root/'outputs/TCR_boundary_independent_audit_v1.json';save(op,report);save(root/'provenance/terminal_v1.json',dict(report,output_SHA256=sha(op)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
