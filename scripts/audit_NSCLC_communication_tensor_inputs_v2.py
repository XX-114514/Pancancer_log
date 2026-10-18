"""Reconstruct every blind partition, author-type gate and native LR position."""
import hashlib
import inspect
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2)
    src=root/'NSCLC_communication_tensor_inputs_20261019_v1'
    out=root/'NSCLC_communication_tensor_inputs_audit_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((src/'private_terminal.json').read_text())
    assert t['source_sha256']=='a6beeadcc7f750c899e2bf237b67248888568de16b8973064017ca26576720fd'
    for name,record in t['files'].items():assert sha(src/name)==record['sha256']
    assert sha(root/'NSCLC_original_author_lineage_audit_20261019_v1/private_terminal.json')=='ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d'
    assert sha(root/'NSCLC_immune_RNA_MCP_clinical_scope_audit_20261019_v1/private_terminal.json')=='78d83326e28986f614de71d3fe3883424e7a3fc850642673f119c9a9af0d78d4'
    read=lambda p:pd.read_csv(p,sep='\t',keep_default_na=False)
    u=read(root/'NSCLC_full_native_gene_units_20261016_v1/private_all_native_units_and_quality.tsv')
    c=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID')
    fp=read(root/'NSCLC_original_author_lineage_qualification_20261019_v1/private_all243_nonoutcome_fingerprint_diagnostic.tsv').set_index('sampleID').private_nonoutcome_fingerprint
    b=read(src/'private_all243_blind_record_partitions.tsv').set_index('sampleID')
    assert list(b.columns)==['profile_block','cancer_type','agent_identity','native_immune_cells','partition','exclusion_reason','reviewed_independent_patient']
    assert len(b)==len(c)==243 and set(b.index)==set(c.index) and not b.reviewed_independent_patient.any()
    c=c.reindex(b.index);n=u.groupby('sampleID').n_cells.sum().reindex(b.index)
    agents=c['anti-PD1_therapy'].apply(lambda x:re.sub(r'\(C\d+(?:-C\d+)?\)','',x).strip())
    base=c.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE')&c.nonoutcome_controls_complete.str.lower().eq('true')&n.ge(100)
    counts=agents[base].value_counts();eligible=base&agents.map(counts).fillna(0).ge(10)
    expected=pd.Series('EXCLUDED',index=b.index)
    cc=pd.DataFrame({'cancer':c.cancer_type,'agent':agents,'profile':fp.reindex(b.index)})
    for _,g in cc[eligible].groupby(['cancer','agent']):
        blocks=list(g.profile.unique());blocks.sort(key=lambda x:hashlib.sha256(('2026101910|'+x).encode()).hexdigest())
        count=len(blocks);nd=count if count==1 else max(1,min(count-1,math.floor(.6*count+.5)))
        train=set(blocks[:nd]);expected.loc[g.index]=np.where(g.profile.isin(train),'DISCOVERY','HELDOUT')
    np.testing.assert_array_equal(b.partition,expected);np.testing.assert_array_equal(b.profile_block,cc.profile)
    np.testing.assert_array_equal(b.agent_identity,agents);np.testing.assert_array_equal(b.cancer_type,c.cancer_type);np.testing.assert_array_equal(b.native_immune_cells,n)
    reasons=[]
    for i in b.index:
        reasons.append('NOT_STRICT_EXPLICIT_ICI' if c.at[i,'strict_treatment_role']!='EXPLICIT_ICI_AUTHOR_CASE' else 'MISSING_FIXED_NONOUTCOME_CONTROLS' if c.at[i,'nonoutcome_controls_complete'].lower()!='true' else 'FEWER_THAN100_IMMUNE_CELLS' if n[i]<100 else 'AGENT_IDENTITY_FEWER_THAN10' if not eligible[i] else 'ELIGIBLE')
    np.testing.assert_array_equal(b.exclusion_reason,reasons)
    assert b[b.partition.ne('EXCLUDED')].groupby('profile_block').partition.nunique().eq(1).all()
    train=b.index[b.partition.eq('DISCOVERY')];test=b.index[b.partition.eq('HELDOUT')]
    s=read(src/'private_all51_author_type_coverage.tsv').set_index('original_author_type')
    assert len(s)==51 and set(s.index)==set(u.sub_cell_type)
    for label,row in s.iterrows():
        g=u[u.sub_cell_type.eq(label)];by=g.set_index('sampleID').n_cells
        nt=int(by.reindex(train,fill_value=0).ge(20).sum());nh=int(by.reindex(test,fill_value=0).ge(20).sum());ok=nt>=math.ceil(len(train)/3)
        assert row.major==g.major_cell_type.iloc[0] and g.major_cell_type.nunique()==1
        assert [row.all_native_cells,row.discovery_contexts,row.discovery_with20cells,row.heldout_contexts,row.heldout_with20cells]==[g.n_cells.sum(),len(train),nt,len(test),nh]
        assert row.eligible==ok and row.status==('DISCOVERY_COVERAGE_ELIGIBLE' if ok else 'DISCOVERY_COVERAGE_BELOW_ONE_THIRD')
    import liana
    resource_path=Path(inspect.getfile(liana.resource.select_resource)).parent/'omni_resource.csv'
    original=pd.read_csv(resource_path);full=original[original.resource.eq('consensus')][['source_genesymbol','target_genesymbol']].rename(columns={'source_genesymbol':'ligand','target_genesymbol':'receptor'}).reset_index(drop=True)
    exported=read(src/'private_complete_original_consensus_resource.tsv')
    pd.testing.assert_frame_equal(full,exported,check_dtype=False)
    genes=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().tolist();lookup={g:i for i,g in enumerate(genes)};assert len(lookup)==31831
    raw=read(src/'private_all_original_LR_resource_positions.tsv')
    unique=json.loads((src/'private_all_unique_LR_complex_positions.json').read_text());byid={r['canonical_resource_LR_id']:r for r in unique};assert len(byid)==len(unique)
    def parts(token):
        token=str(token).strip().upper()
        if token in lookup:return [token]
        return sorted(set(token.split('_'))) if token and token!='NAN' else []
    requested=set();seen=set()
    assert len(raw)==len(full)
    for i,(ligand,receptor) in enumerate(full.itertuples(index=False,name=None)):
        L=parts(ligand);R=parts(receptor);key=json.dumps([L,R],separators=(',',':'));rid='LR_'+hashlib.sha256(key.encode()).hexdigest()[:16]
        missing=sorted(x for x in set(L+R) if x not in lookup);contested=bool(set(L)&{'C1QA','C1QB','C1QC'} and 'CD93' in R);ok=bool(L and R and not missing and not contested)
        status='PRIMARY_NATIVE_COMPLEX_AVAILABLE' if ok else 'DESCRIPTIVE_CONTESTED_C1Q_CD93_BINDING' if contested else 'MISSING_NATIVE_COMPLEX_COMPONENT' if missing else 'EMPTY_RESOURCE_COMPONENT'
        r=raw.iloc[i];assert [r.resource_row,r.ligand_complex_original,r.receptor_complex_original,r.canonical_resource_LR_id,r.metadata_primary_eligible,r.status]==[i,ligand,receptor,rid,ok,status]
        d=byid[rid];assert d==dict(canonical_resource_LR_id=rid,canonical_component_key=key,ligand_components=L,receptor_components=R,native_ligand_columns=[lookup.get(g,-1) for g in L],native_receptor_columns=[lookup.get(g,-1) for g in R],missing_native_components=missing,metadata_primary_eligible=ok,binding_contested=contested,status=status)
        seen.add(rid)
        if ok:requested.update(L+R)
    assert seen==set(byid)
    axis=read(src/'private_complete_primary_LR_gene_native_axis.tsv');wanted=sorted(requested,key=lookup.get)
    assert axis.gene.tolist()==wanted and axis.native_column.tolist()==[lookup[g] for g in wanted]
    coverage=read(src/'private_all_source_blind_partition_coverage.tsv').set_index(['partition','cancer_type','agent_identity'])
    assert len(coverage)==b.groupby(['partition','cancer_type','agent_identity']).ngroups
    for key,g in b.groupby(['partition','cancer_type','agent_identity']):
        r=coverage.loc[key];assert r.source=='GSE243013' and r.author_records==len(g) and r.profile_blocks==g.profile_block.nunique() and r.reviewed_independent_patients==0
    reg=read(src/'private_all20_prospective_factor_clinical_slots.tsv');assert len(reg)==20 and not reg.model_id.duplicated().any()
    for r,(k,f) in zip(reg.itertuples(),[(k,f) for k in range(2,7) for f in range(1,k+1)]):
        assert [r.model_id,r.rank,r.factor,r.training_restarts,r.clinical_partition,r.family_denominator,r.outcome,r.status]==['TC'+str(k).zfill(2)+'F'+str(f).zfill(2),k,f,3,'HELDOUT',20,'author_pCR_or_MPR_1_vs_nonMPR_0','FUTURE_FACTOR_SLOT_NO_SCORES_OR_EFFECT']
    p=json.loads((src/'private_resource_and_scope_provenance.json').read_text())
    assert p['bundled_full_LIANA_resource_sha256']==sha(resource_path) and p['source_gene_statistics_not_yet_read'] and p['author_type_count_corrected_to_raw51']
    assert p['maximum_discovery_tensor_float32_bytes_before_expression_support_filter']==len(train)*sum(x['metadata_primary_eligible'] for x in unique)*int(s.eligible.sum())**2*4
    receipt={'status':'ACCEPTED_ALL243_BLIND_PARTITIONS_ALL51_TYPE_GATES_COMPLETE_RESOURCE_NATIVE_LR_AXIS_AND20_SLOTS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'discovery_records':len(train),'heldout_records':len(test),'eligible_types':int(s.eligible.sum()),'primary_LR_gene_columns':len(axis),'new_expression_or_effects':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
