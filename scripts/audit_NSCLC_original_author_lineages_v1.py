"""Independent raw-metadata replay for every lineage row, pool and record."""
import hashlib
import json
import os
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import FIELDS
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_original_author_lineage_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'NSCLC_original_author_lineage_qualification_20261019_v1';t=wait(root,out,'qualify_NSCLC_original_author_lineages_v1',src/'private_terminal.json',{'COMPLETE_ALL_NATIVE_AUTHOR_LINEAGE_AND_RECORD_QUALIFICATION_AUDIT_PENDING'})
    for name,v in t['files'].items():assert sha(src/name)==v['sha256']
    native=root/'NSCLC_full_native_gene_units_20261016_v1';unit=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).sort_values('unit_index');oldmap=np.load(native/'private_cell_unit_index.npy',mmap_mode='r');RNA=np.load(native/'private_full_RNA.npy',mmap_mode='r')
    unitlookup={(r.sampleID,r.sub_cell_type):(r.unit_index,r.major_cell_type,r.n_cells) for r in unit.itertuples(index=False)};assert len(unitlookup)==11701
    old=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1';selected=pd.read_csv(old/'private_selected_native_cells.tsv.gz',sep='\t');selected_lookup={int(r.native_row):(r.sampleID,r.view) for r in selected.itertuples(index=False)};assert len(selected_lookup)==len(selected)
    legacy=pd.read_csv(old/'private_all729_material_view_capture.tsv',sep='\t').set_index(['sampleID','view']);views=['R_FOXP3_capture','E_effector_capture','Myeloid']
    clinical={};unit_counts=Counter();categoryrows=defaultdict(list);purity=Counter();all_labels=defaultdict(set);row=0
    cols=list(dict.fromkeys(['cellID','sampleID','sub_cell_type','major_cell_type','total_counts']+FIELDS));barcode=pd.read_csv(root/'GSE243013_barcodes.csv.gz',dtype=str,keep_default_na=False,chunksize=125000)
    for chunk,bcs in zip(pd.read_csv(root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',usecols=cols,dtype=str,keep_default_na=False,chunksize=125000),barcode):
        assert len(chunk)==len(bcs) and np.array_equal(chunk.cellID,bcs.barcode)
        order=chunk[['sampleID','sub_cell_type','major_cell_type','total_counts']]
        for sid,sub,major,total in order.itertuples(index=False,name=None):
            ui,expected_major,_=unitlookup[(sid,sub)];assert expected_major==major and ui==int(oldmap[row]) and int(float(total))==int(RNA[row]);unit_counts[(sid,sub)]+=1
            side='R_FOXP3_capture' if major=='T/NK cell' and sub.startswith('CD4T_Treg_') else 'E_effector_capture' if major=='T/NK cell' and sub.startswith('CD8T_') else 'Myeloid' if major=='Myeloid cell' else None
            if side is not None:categoryrows[(sid,side)].append(row);all_labels[side].add(sub)
            if row in selected_lookup:
                ss,v=selected_lookup[row];assert ss==sid;purity[(sid,v,sub,side==v)]+=1
            row+=1
        for record in chunk[FIELDS].drop_duplicates().itertuples(index=False,name=None):
            if record[0] in clinical:assert clinical[record[0]]==record
            else:clinical[record[0]]=record
    assert row==1254749 and len(clinical)==243
    for key,(_,_,n) in unitlookup.items():assert unit_counts[key]==n
    c=pd.read_csv(src/'private_all243_original_author_clinical_records.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False);expected=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False).loc[c.index];pd.testing.assert_frame_equal(c,expected)
    for sid in c.index:assert tuple(c.loc[sid,FIELDS])==clinical[sid]
    pools=pd.read_csv(src/'private_all729_author_lineage_pool_positions.tsv',sep='\t');assert len(pools)==729 and pools.pool_index.tolist()==list(range(729))
    expected_members=[];expected_selected=[];samples=sorted(clinical)
    for p in pools.itertuples(index=False):
        si,vi=divmod(p.pool_index,3);sid,view=samples[si],views[vi];assert p.sampleID==sid and p.view==view and p.native_material_index==si and p.view_index==vi
        rows=categoryrows[(sid,view)];n=len(rows);oldn=int(legacy.loc[(sid,view),'full_view_cells']);eligible=n>=20 and oldn>=20
        assert p.all_author_cells==n and p.original_capture_cells==oldn and p.eligible==eligible and p.all_native_material_cells==int(legacy.loc[(sid,view),'all_native_material_cells'])
        reason='AUTHOR_AND_ORIGINAL_CAPTURE_ELIGIBLE' if eligible else 'AUTHOR_FEWER_THAN20' if n<20 else 'ORIGINAL_CAPTURE_FEWER_THAN20';assert p.status==reason and p.patient_independence_status=='AUTHOR_RECORD_ONLY_NO_REVIEWED_PATIENT_ALIAS'
        assert p.selected_cells==(min(100,n) if eligible else 0)
        if eligible:
            order=sorted((hashlib.sha256(('2026101802|'+view+'|'+str(rr)).encode()).digest(),rr) for rr in rows)
            chosen={rr for _,rr in order[:100]}
            expected_members.extend((rr,p.pool_index,rr in chosen) for rr in rows)
    expected_members.sort();a=dict(np.load(src/'private_all_author_lineage_native_memberships.npz'));assert np.array_equal(a['native_row'],[x[0] for x in expected_members]) and np.array_equal(a['pool_index'],[x[1] for x in expected_members]) and np.array_equal(a['selected'],[x[2] for x in expected_members])
    table=pd.read_csv(src/'private_all_cached_marker_capture_author_reference_counts.tsv',sep='\t');actual={(r.sampleID,r.view,r.original_author_subtype,r.agrees_with_author_lineage):int(r.selected_cells) for r in table.itertuples(index=False)};assert len(actual)==len(table) and actual==purity
    labels=json.loads((src/'private_original_reference_labels_and_identity_limits.json').read_text());assert labels['labels']=={v:sorted(s) for v,s in all_labels.items()} and labels['explicit_patient_columns']==[] and labels['no_alias_adopted'] and not labels['formal_state_changed']
    nonoutcome=[x for x in FIELDS if x not in ['sampleID','pathological_response','pathological_response_rate','radiological_response']];finger={sid:hashlib.sha256(json.dumps([c.loc[sid,k] for k in nonoutcome],ensure_ascii=False,separators=(',',':')).encode()).hexdigest() for sid in c.index};counts=Counter(finger.values());f=pd.read_csv(src/'private_all243_nonoutcome_fingerprint_diagnostic.tsv',sep='\t')
    assert len(f)==243
    for r in f.itertuples(index=False):assert r.private_nonoutcome_fingerprint==finger[r.sampleID] and r.same_nonoutcome_records==counts[finger[r.sampleID]] and r.not_a_patient_alias
    assert labels['nonoutcome_fingerprint_unique']==len(counts) and labels['records_in_duplicate_nonoutcome_profiles']==sum(v for v in counts.values() if v>1) and labels['profile_groups_above_one']==sum(v>1 for v in counts.values())
    coverage=pd.read_csv(src/'private_all_source_role_author_lineage_pair_coverage.tsv',sep='\t');expected_keys={(role,cancer,va,vb) for role,cancer in set(zip(c.strict_treatment_role,c.cancer_type)) for va,vb in [(views[0],views[1]),(views[0],views[2]),(views[1],views[2])]};actual_keys=set()
    for r in coverage.itertuples(index=False):
        key=(r.treatment_role,r.cancer_type,r.view_A,r.view_B);assert key not in actual_keys;actual_keys.add(key)
        cohort=c[(c.strict_treatment_role==r.treatment_role)&(c.cancer_type==r.cancer_type)];eligible=[]
        for sid in cohort.index:
            if all(len(categoryrows[(sid,v)])>=20 and int(legacy.loc[(sid,v),'full_view_cells'])>=20 for v in [r.view_A,r.view_B]):eligible.append(sid)
        use=c.loc[eligible];good=use.nonoutcome_controls_complete.str.lower().eq('true');path=use.author_pathology_evaluable.str.lower().eq('true')
        assert r.all_author_records==len(cohort) and r.both_author_and_legacy_captures==len(eligible) and r.complete_old_nonoutcome_controls==int(good.sum()) and r.complete_controls_and_pathology==int((good&path).sum()) and r.reviewed_independent_patient_count==0
    assert actual_keys==expected_keys
    receipt={'status':'ACCEPTED_ALL1254749_AUTHOR_ROWS729_POOLS_ALL_MEMBERSHIPS_HASHES243_RECORDS_AND_CAPTURE_REFERENCE_COUNTS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'native_metadata_rows':row,'author_records':243,'all_pool_positions':729,'all_eligible_author_cells':len(expected_members),'selected_cells':int(a['selected'].sum()),'new_patient_aliases_or_expression_or_effects':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
