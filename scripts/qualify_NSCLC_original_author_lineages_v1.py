"""Full native author-lineage/record qualification before any new RNA assay."""
import hashlib
import json
import os
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import FIELDS,COMMON


VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
TREG={'CD4T_Treg_CCR8','CD4T_Treg_FOXP3','CD4T_Treg_MKI67'}


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    assertions={'NSCLC_immune_RNA_MCP_native_cells_audit_20261019_v1':'de1f3ede3b2a62559c368f26a30507e1dcb676a0ebb343a4d82ec7b69780e210','NSCLC_immune_RNA_MCP_clinical_scope_audit_20261019_v1':'78d83326e28986f614de71d3fe3883424e7a3fc850642673f119c9a9af0d78d4','NSCLC_full_native_gene_unit_audit_20261016_v1':'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'}
    for folder,digest in assertions.items():assert sha(root/folder/'private_terminal.json')==digest
    source=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz';assert sha(source)=='24155466454df53d01eb6cbbfcd56abeb45adb68f74e721af22a9dce41040b2f'
    out=root/'NSCLC_original_author_lineage_qualification_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    header=pd.read_csv(source,nrows=0).columns.tolist()
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'prerequisite_acceptance_sha256':assertions,'original_metadata_sha256':sha(source),'scope':'Every1254749 original native cells/all243 author records/all11701 original units. No expression read, programme association, learnedgene selection or clinical endpoint change. Original labels not adopted as new global states.','lineages':{'R':'All3 exact original CD4T_Treg labels','E':'Every original CD8T_ prefix category, including author-labelled MAIT/NK-like categories. This is an author CD8 union, not proof of conventional-only CD8 equivalence across sources.','My':'All original Myeloid cell major labels'},'coverage':'All243x3positions retained. Prospective hash100 from allauthor lineage only if both author>=20 and corresponding original frozen marker/major capture>=20. Preserve old capture fraction/logN meaning, unknown source offset, original allimmune denominator and training-domain caveat.','identity':'Exact original row/metadata/barcode and original unit keys. Unique sampleID alone is an author record, not reviewed patient identity. No inferred alias or grouping by clinical fingerprint. Identical nonoutcome clinical fingerprints flagged as sensitivity opportunities only.','source_context':'Official GSE243013/BioProject source protocol postneoadjuvant surgical material;234study count versus243deposition distinction unbridged. Current cell metadata has no separate explicit patient key. No perrecord phase or independent patient status promotion.','capture_comparison':'Existing hash-selected original marker captures only: complete author category counts for every cached selected cell. Not whole-capture purity or formal doublet detection.','future_scope':'Full two frozen MCPs and all three available lineages should be considered in new context; all50 reference qualification would need exclude predictor union including Myeloid. No new900-family registration/effects here.'}
    (out/'plan_before_author_lineage_or_new_assay.json').write_text(json.dumps(plan,indent=2)+'\n')
    native=root/'NSCLC_full_native_gene_units_20261016_v1';unit=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).sort_values('unit_index');assert unit.unit_index.tolist()==list(range(11701))
    original_units=np.load(native/'private_cell_unit_index.npy',mmap_mode='r');RNA=np.load(native/'private_full_RNA.npy',mmap_mode='r');assert len(original_units)==1254749
    unitkeys=pd.MultiIndex.from_frame(unit[['sampleID','sub_cell_type']]);major=unit.major_cell_type.to_numpy(str);sub=unit.sub_cell_type.to_numpy(str);usample=unit.sampleID.to_numpy(str)
    native_view=np.full(len(original_units),-1,np.int8);clinical={};offset=0;unitcounts=np.zeros(len(unit),np.int64)
    barcode_chunks=pd.read_csv(root/'GSE243013_barcodes.csv.gz',dtype=str,keep_default_na=False,chunksize=200000)
    fields=list(dict.fromkeys(['cellID','sampleID','sub_cell_type','major_cell_type','total_counts']+FIELDS))
    for chunk,bcs in zip(pd.read_csv(source,usecols=fields,dtype=str,keep_default_na=False,chunksize=200000),barcode_chunks):
        n=len(chunk);assert n==len(bcs) and np.array_equal(chunk.cellID,bcs.barcode)
        index=unitkeys.get_indexer(pd.MultiIndex.from_frame(chunk[['sampleID','sub_cell_type']]))
        assert (index>=0).all() and np.array_equal(index,original_units[offset:offset+n])
        assert np.array_equal(major[index],chunk.major_cell_type) and np.array_equal(usample[index],chunk.sampleID)
        assert np.array_equal(pd.to_numeric(chunk.total_counts).to_numpy(np.int64),RNA[offset:offset+n])
        unitcounts+=np.bincount(index,minlength=len(unit))
        m=chunk.major_cell_type.to_numpy(str);s=chunk.sub_cell_type.to_numpy(str);v=np.full(n,-1,np.int8);v[(m=='T/NK cell')&np.isin(s,list(TREG))]=0;v[(m=='T/NK cell')&np.char.startswith(s,'CD8T_')]=1;v[m=='Myeloid cell']=2;native_view[offset:offset+n]=v
        for record in chunk[FIELDS].drop_duplicates().itertuples(index=False,name=None):
            sid=record[0]
            if sid in clinical:assert clinical[sid]==record
            else:clinical[sid]=record
        offset+=n
    assert offset==1254749 and len(clinical)==243 and np.array_equal(unitcounts,unit.n_cells)
    assert set(sub[np.char.startswith(sub,'CD4T_Treg_')])==TREG
    labels={view:sorted(set(sub[(unit.major_cell_type.eq('T/NK cell')&unit.sub_cell_type.isin(TREG)) if j==0 else (unit.major_cell_type.eq('T/NK cell')&unit.sub_cell_type.str.startswith('CD8T_')) if j==1 else unit.major_cell_type.eq('Myeloid cell')])) for j,view in enumerate(VIEWS)}
    assert len(labels[VIEWS[1]])==10
    original=pd.DataFrame(list(clinical.values()),columns=FIELDS).sort_values('sampleID').set_index('sampleID',drop=False)
    c=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False).loc[original.index]
    for field in FIELDS:assert np.array_equal(original[field].astype(str),c[field].astype(str))
    assert not c.reviewed_independent_patient.astype(str).str.lower().eq('true').any()
    samples=original.index.tolist();samplemap={x:i for i,x in enumerate(samples)};sample_index=np.array([samplemap[x] for x in usample],np.int16)[original_units]
    old=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1';legacy=pd.read_csv(old/'private_all729_material_view_capture.tsv',sep='\t').set_index(['sampleID','view']);selected=pd.read_csv(old/'private_selected_native_cells.tsv.gz',sep='\t');assert not selected.native_row.duplicated().any()
    purity=[];counts=[]
    for r in selected.itertuples(index=False):
        row=int(r.native_row);ui=int(original_units[row]);assert usample[ui]==r.sampleID
        purity.append({'sampleID':r.sampleID,'view':r.view,'original_author_subtype':sub[ui],'agrees_with_author_lineage':bool(native_view[row]==VIEWS.index(r.view)),'selected_cells':1})
    purity=pd.DataFrame(purity).groupby(['sampleID','view','original_author_subtype','agrees_with_author_lineage']).selected_cells.sum().reset_index()
    purity.to_csv(out/'private_all_cached_marker_capture_author_reference_counts.tsv',sep='\t',index=False)
    group=sample_index.astype(np.int64)*3+np.maximum(native_view,0);valid=native_view>=0;fullcount=np.bincount(group[valid],minlength=243*3).reshape(243,3)
    memberships=[];poolrecords=[];pid=0
    for i,sid in enumerate(samples):
        for j,view in enumerate(VIEWS):
            oldcount=int(legacy.loc[(sid,view),'full_view_cells']);cells=np.flatnonzero((sample_index==i)&(native_view==j));n=int(fullcount[i,j]);assert n==len(cells)
            eligible=n>=20 and oldcount>=20
            status='AUTHOR_AND_ORIGINAL_CAPTURE_ELIGIBLE' if eligible else 'AUTHOR_FEWER_THAN20' if n<20 else 'ORIGINAL_CAPTURE_FEWER_THAN20'
            poolrecords.append({'pool_index':pid,'sampleID':sid,'view':view,'view_index':j,'native_material_index':i,'all_author_cells':n,'original_capture_cells':oldcount,'all_native_material_cells':int(legacy.loc[(sid,view),'all_native_material_cells']),'eligible':eligible,'status':status,'selected_cells':min(100,n) if eligible else 0,'patient_independence_status':'AUTHOR_RECORD_ONLY_NO_REVIEWED_PATIENT_ALIAS'})
            if eligible:
                chosen=set(sorted(cells,key=lambda rr:hashlib.sha256(('2026101802|'+view+'|'+str(int(rr))).encode()).hexdigest())[:100])
                memberships.extend((int(rr),pid,rr in chosen) for rr in cells)
            pid+=1
    pool=pd.DataFrame(poolrecords);pool.to_csv(out/'private_all729_author_lineage_pool_positions.tsv',sep='\t',index=False)
    memberships.sort();rows=np.array([x[0] for x in memberships],np.int32);assert len(rows)==len(set(rows))
    np.savez_compressed(out/'private_all_author_lineage_native_memberships.npz',native_row=rows,pool_index=np.array([x[1] for x in memberships],np.int16),selected=np.array([x[2] for x in memberships],bool))
    c.to_csv(out/'private_all243_original_author_clinical_records.tsv',sep='\t',index=False)
    # Fingerprints describe indistinguishable deposited nonoutcome fields only.
    nonoutcome=[x for x in FIELDS if x not in ['sampleID','pathological_response','pathological_response_rate','radiological_response']]
    fingerprint=c[nonoutcome].astype(str).apply(lambda rr:hashlib.sha256(json.dumps(rr.tolist(),ensure_ascii=False,separators=(',',':')).encode()).hexdigest(),axis=1)
    sizes=fingerprint.value_counts();pd.DataFrame({'sampleID':c.sampleID,'private_nonoutcome_fingerprint':fingerprint,'same_nonoutcome_records':fingerprint.map(sizes),'not_a_patient_alias':True}).to_csv(out/'private_all243_nonoutcome_fingerprint_diagnostic.tsv',sep='\t',index=False)
    coverage=[]
    for (role,cancer),g in c.groupby(['strict_treatment_role','cancer_type'],sort=True):
        ids=set(g.sampleID)
        for va,vb in [(VIEWS[0],VIEWS[1]),(VIEWS[0],VIEWS[2]),(VIEWS[1],VIEWS[2])]:
            A=pool[pool.view.eq(va)&pool.eligible].set_index('sampleID');B=pool[pool.view.eq(vb)&pool.eligible].set_index('sampleID');common=sorted(ids&set(A.index)&set(B.index));usable=c.loc[common]
            complete=usable.nonoutcome_controls_complete.astype(str).str.lower().eq('true');pathology=usable.author_pathology_evaluable.astype(str).str.lower().eq('true')
            coverage.append({'source':'GSE243013','treatment_role':role,'cancer_type':cancer,'view_A':va,'view_B':vb,'all_author_records':len(g),'both_author_and_legacy_captures':len(common),'complete_old_nonoutcome_controls':int(complete.sum()),'complete_controls_and_pathology':int((complete&pathology).sum()),'reviewed_independent_patient_count':0,'source_context':'post_neoadjuvant_surgical; individual_phase_not_promoted'})
    pd.DataFrame(coverage).to_csv(out/'private_all_source_role_author_lineage_pair_coverage.tsv',sep='\t',index=False)
    (out/'private_original_reference_labels_and_identity_limits.json').write_text(json.dumps({'labels':labels,'original_metadata_columns':header,'explicit_patient_columns':[x for x in header if 'patient' in x.lower() or 'donor' in x.lower()],'all243_original_clinical_records_identical':True,'nonoutcome_fingerprint_unique':int(len(sizes)),'records_in_duplicate_nonoutcome_profiles':int(sizes[sizes>1].sum()),'profile_groups_above_one':int((sizes>1).sum()),'no_alias_adopted':True,'post_neoadjuvant_context_source_only':True,'author_CD8_union_includes_MAIT':True,'formal_state_changed':False},indent=2)+'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt={'status':'COMPLETE_ALL_NATIVE_AUTHOR_LINEAGE_AND_RECORD_QUALIFICATION_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'original_native_cells':offset,'original_units':len(unit),'author_records':243,'registered_pool_positions':729,'eligible_pools':int(pool.eligible.sum()),'full_eligible_author_cells':len(rows),'selected_cells':sum(x[2] for x in memberships),'new_expression_or_associations':False,'patient_aliases_adopted':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
