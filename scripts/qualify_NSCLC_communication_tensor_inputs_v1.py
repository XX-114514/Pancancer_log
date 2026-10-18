"""Outcome-blind blocks, full author labels and complete native LR vocabulary."""
import hashlib
import inspect
import json
import math
import os
import re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    assert sha(root/'NSCLC_original_author_lineage_audit_20261019_v1/private_terminal.json')=='ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d'
    assert sha(root/'NSCLC_immune_RNA_MCP_clinical_scope_audit_20261019_v1/private_terminal.json')=='78d83326e28986f614de71d3fe3883424e7a3fc850642673f119c9a9af0d78d4'
    out=root/'NSCLC_communication_tensor_inputs_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'scope':'All original51author categories/all243records/alluniqueLIANAconsensus LRcomplexes. Exact direct native symbol components only; no endpoint chosen by programme/response effects. Source immune-only, not stromal/epithelial or original9CM validation.','partition':'Strict original explicitICI role, complete original8nonoutcome controls,>=100nativeimmune cells; drug-identity groups>=10source records. Profile block uses all original nonoutcome clinical fields only, no patient alias assertion. Same profile stays in onefold. Per cancer x canonicaldrug stratum hash2026101910 profile sorting; nearest60percent discovery blocks constrained1..n-1 when n>=2,singleton discovery. No outcome balancing, seed search or rerandomization.','states':'All51registered. Tensor state requires >=20cells in at leastceil(discoveryN/3) discovery contexts; missingstate remains mask0, actual measured expressionzero stays observed. Heldout not used for state/resource/support selection.','native_LR':'LIANA1.7.1 bundledconsensus, unique molecular-component sets deduplicated preserving orientation. Every complex component required on native31831axis; missing/ambiguous neverzero. C1Qligand toCD93receptor retained descriptive and excluded primarybinding vocabulary per prior project evidence rule.','future_scoring':'All-cell original-unit meanlogCP10K, complexminimum across subunits, exprprod and log1p transform. LRsupport>=10percent componentdetection in some eligible sender AND receiver state in at leastceil(discoveryN/3) discovery contexts. This is a global discovery-only feature filter; all percontext measurable values retained even below10percent. Only unit<20/geneunmeasurable is missing, not an interaction missing by significance.','future_model':'Official Tensor-cell2cell/PARAFAC, allranks2..6 and3fixedrandomstarts/rank; best observed discovery reconstruction only, all20factor clinical slots retained. Missingmask/knownzero and fixedfactor heldout projection must pass synthetic gate. No clinical scores/truefactor fit executed here. Detailed actualAPI/optimization/resource contract must be frozen before fit.','limits':'Authorrecords not independently reviewed patients, wholeproject source/outcomes previously explored. Source postneoadjuvant surgery association, not baselineprediction. No communication strength/binding/causal claim from exprprod or latent factors. OldMCP/400functional clinical families remain closed.'}
    (out/'plan_before_context_split_resource_or_tensor.json').write_text(json.dumps(plan,indent=2)+'\n')
    native=root/'NSCLC_full_native_gene_units_20261016_v1';units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False);assert len(units)==11701 and units.sub_cell_type.nunique()==51
    c=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False).sort_index();assert len(c)==243
    fingerprint=pd.read_csv(root/'NSCLC_original_author_lineage_qualification_20261019_v1/private_all243_nonoutcome_fingerprint_diagnostic.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID').private_nonoutcome_fingerprint
    total=units.groupby('sampleID').n_cells.sum();c['native_immune_cells']=total.reindex(c.index);c['agent_identity']=c['anti-PD1_therapy'].map(lambda s:re.sub(r'\(C\d+(?:-C\d+)?\)','',s).strip());c['profile_block']=fingerprint.reindex(c.index)
    base=c.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE')&c.nonoutcome_controls_complete.str.lower().eq('true')&c.native_immune_cells.ge(100);agentcounts=c.loc[base,'agent_identity'].value_counts();eligible=base&c.agent_identity.isin(agentcounts[agentcounts>=10].index)
    c['partition']='EXCLUDED';c['exclusion_reason']=np.where(~c.strict_treatment_role.eq('EXPLICIT_ICI_AUTHOR_CASE'),'NOT_STRICT_EXPLICIT_ICI',np.where(~c.nonoutcome_controls_complete.str.lower().eq('true'),'MISSING_FIXED_NONOUTCOME_CONTROLS',np.where(c.native_immune_cells.lt(100),'FEWER_THAN100_IMMUNE_CELLS',np.where(~eligible,'AGENT_IDENTITY_FEWER_THAN10','ELIGIBLE'))))
    for _,g in c[eligible].groupby(['cancer_type','agent_identity'],sort=True):
        blocks=sorted(set(g.profile_block),key=lambda b:hashlib.sha256(('2026101910|'+b).encode()).hexdigest());n=len(blocks);nd=n if n==1 else min(n-1,max(1,int(math.floor(.6*n+.5))));discovery=set(blocks[:nd]);c.loc[g.index,'partition']=['DISCOVERY' if b in discovery else 'HELDOUT' for b in g.profile_block]
    assert c[c.partition.ne('EXCLUDED')].groupby('profile_block').partition.nunique().eq(1).all()
    blind=c[['sampleID','profile_block','cancer_type','agent_identity','native_immune_cells','partition','exclusion_reason']].copy();blind['reviewed_independent_patient']=False;blind.to_csv(out/'private_all243_blind_record_partitions.tsv',sep='\t',index=False)
    discovery=c.index[c.partition.eq('DISCOVERY')].tolist();heldout=c.index[c.partition.eq('HELDOUT')].tolist();assert len(discovery)>0 and len(heldout)>0
    types=sorted(units.sub_cell_type.unique());state=[]
    for label in types:
        g=units[units.sub_cell_type.eq(label)];assert g.major_cell_type.nunique()==1;counts=g.set_index('sampleID').n_cells;nt=int(counts.reindex(discovery,fill_value=0).ge(20).sum());nh=int(counts.reindex(heldout,fill_value=0).ge(20).sum());ok=nt>=math.ceil(len(discovery)/3)
        state.append({'original_author_type':label,'major':g.major_cell_type.iloc[0],'all_native_cells':int(g.n_cells.sum()),'discovery_contexts':len(discovery),'discovery_with20cells':nt,'heldout_contexts':len(heldout),'heldout_with20cells':nh,'eligible':ok,'status':'DISCOVERY_COVERAGE_ELIGIBLE' if ok else 'DISCOVERY_COVERAGE_BELOW_ONE_THIRD'})
    pd.DataFrame(state).to_csv(out/'private_all51_author_type_coverage.tsv',sep='\t',index=False)
    import liana
    selected=liana.resource.select_resource('consensus').copy();assert set(selected.columns)=={'ligand','receptor'};selected.to_csv(out/'private_complete_original_consensus_resource.tsv',sep='\t',index=False)
    source_csv=Path(inspect.getfile(liana.resource.select_resource)).parent/'omni_resource.csv';assert source_csv.exists()
    genes=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().tolist();lookup={g:i for i,g in enumerate(genes)};assert len(lookup)==31831
    def parse(token):
        token=str(token).strip().upper();return [token] if token in lookup else sorted(set(token.split('_'))) if token and token!='NAN' else []
    raw=[];unique={};requested=set()
    for i,(ligand,receptor) in enumerate(selected.itertuples(index=False,name=None)):
        L,R=parse(ligand),parse(receptor);key=json.dumps([L,R],separators=(',',':'));rid='LR_'+hashlib.sha256(key.encode()).hexdigest()[:16];missing=sorted(g for g in set(L+R) if g not in lookup);contested=bool('CD93' in R and {'C1QA','C1QB','C1QC'}&set(L));eligible=bool(L and R and not missing and not contested)
        status='PRIMARY_NATIVE_COMPLEX_AVAILABLE' if eligible else 'DESCRIPTIVE_CONTESTED_C1Q_CD93_BINDING' if contested else 'MISSING_NATIVE_COMPLEX_COMPONENT' if missing else 'EMPTY_RESOURCE_COMPONENT'
        raw.append({'resource_row':i,'ligand_complex_original':ligand,'receptor_complex_original':receptor,'canonical_resource_LR_id':rid,'metadata_primary_eligible':eligible,'status':status})
        if rid in unique:assert unique[rid]['canonical_component_key']==key
        else:unique[rid]={'canonical_resource_LR_id':rid,'canonical_component_key':key,'ligand_components':L,'receptor_components':R,'native_ligand_columns':[lookup.get(g,-1) for g in L],'native_receptor_columns':[lookup.get(g,-1) for g in R],'missing_native_components':missing,'metadata_primary_eligible':eligible,'binding_contested':contested,'status':status}
        if eligible:requested.update(L+R)
    pd.DataFrame(raw).to_csv(out/'private_all_original_LR_resource_positions.tsv',sep='\t',index=False);details=sorted(unique.values(),key=lambda d:d['canonical_resource_LR_id']);(out/'private_all_unique_LR_complex_positions.json').write_text(json.dumps(details,indent=2)+'\n')
    requested=sorted(requested,key=lookup.get);pd.DataFrame({'gene':requested,'native_column':[lookup[g] for g in requested]}).to_csv(out/'private_complete_primary_LR_gene_native_axis.tsv',sep='\t',index=False)
    coverage=[]
    for (part,cancer,agent),g in c.groupby(['partition','cancer_type','agent_identity'],sort=True):coverage.append({'source':'GSE243013','partition':part,'cancer_type':cancer,'agent_identity':agent,'author_records':len(g),'profile_blocks':g.profile_block.nunique(),'reviewed_independent_patients':0})
    pd.DataFrame(coverage).to_csv(out/'private_all_source_blind_partition_coverage.tsv',sep='\t',index=False)
    reg=[{'model_id':'TC'+str(k).zfill(2)+'F'+str(f).zfill(2),'rank':k,'factor':f,'training_restarts':3,'clinical_partition':'HELDOUT','family_denominator':20,'outcome':'author_pCR_or_MPR_1_vs_nonMPR_0','status':'FUTURE_FACTOR_SLOT_NO_SCORES_OR_EFFECT'} for k in range(2,7) for f in range(1,k+1)];pd.DataFrame(reg).to_csv(out/'private_all20_prospective_factor_clinical_slots.tsv',sep='\t',index=False)
    cap=len(discovery)*sum(d['metadata_primary_eligible'] for d in details)*sum(r['eligible'] for r in state)**2*4
    provenance={'created_UTC':datetime.now(timezone.utc).isoformat(),'bundled_full_LIANA_resource_sha256':sha(source_csv),'liana_version':liana.__version__,'source_genes_sha256':sha(root/'GSE243013_genes.csv.gz'),'native_author_unit_table_sha256':sha(native/'private_all_native_units_and_quality.tsv'),'clinical_scope_sha256':sha(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv'),'profile_diagnostic_sha256':sha(root/'NSCLC_original_author_lineage_qualification_20261019_v1/private_all243_nonoutcome_fingerprint_diagnostic.tsv'),'maximum_discovery_tensor_float32_bytes_before_expression_support_filter':cap,'source_gene_statistics_not_yet_read':True,'author_type_count_corrected_to_raw51':True};(out/'private_resource_and_scope_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt={'status':'COMPLETE_BLIND_RECORD_SPLIT_ALL51_AUTHOR_TYPES_AND_FULL_LR_METADATA_SCOPE_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'original_author_types':51,'original_records':243,'eligible_records':len(discovery)+len(heldout),'discovery_records':len(discovery),'heldout_records':len(heldout),'metadata_eligible_types':sum(r['eligible'] for r in state),'raw_resource_rows':len(raw),'unique_LR_positions':len(details),'metadata_primary_eligible_LRs':sum(d['metadata_primary_eligible'] for d in details),'primary_LR_gene_columns':len(requested),'future_clinical_slots':20,'new_expression_tensor_factor_or_clinical_effects':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
