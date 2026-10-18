"""Outcome-free qualification of all fixed Hallmark complement references."""
import importlib.metadata
import json
import os
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);old=root/'immune_RNA_MCP_functional_description_20261018_v1';assert sha(root/'immune_RNA_MCP_functional_description_audit_20261018_v1/private_terminal.json')=='e96be1d1a2630885557cec7443b931723a8f162f76dbff2b920f92b638f37b42';refs=Path(__file__).parent/'hallmark50_public_reference_v3';manifest=refs/'reference_provenance.json';assert sha(manifest)==json.loads((old/'plan_before_complete_functional_description.json').read_text())['public_reference_manifest_sha256'];records=json.loads(manifest.read_text())['records'];assert len(records)==50
    out=root/'GSE176021_disjoint_function_reference_scope_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    try:version=importlib.metadata.version('pyucell')
    except importlib.metadata.PackageNotFoundError:version='NOT_IN_CURRENT_ENVIRONMENT'
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'reference_manifest_sha256':sha(manifest),'scope':'Every one of the50 already accepted public Hallmark sets. Remove union of all active genes in BOTH frozen programmes and BOTH R/E views to avoid a functional target containing the predictor genes. Full signatures may be descriptive sensitivity, never independent validation. No handpicked pathway or new fitted gene/state.','reference_fidelity_gates':'Keep at least80% of original reference after exclusion; at least80% of remaining genes uniquely source-measurable, at least20 uniquely measured complement genes. All failed positions remain explicit; do not trim programmes or lower gates. Complements are labelled complements, not canonical full signatures.','planned_measurement':'Author-defined lineage scopes, all original measured source-gene universe; per-cell rank-based UCell primary maxrank1500 and1000/2000 fixed sensitivities, mean logCP10K expression sensitivity. Official runtime and independent scoring math must be qualified before expression; no dependency update or new data download here.','prospective_question_to_freeze_after_qualification':'All50references x2programmes x2program-lineages x2response-lineages=400 potential paired-delta association slots; primary same10paired cases, within and cross lineage kept distinct, all eligible slots in one family. This file qualifies references only; no new functional scores or associations.','prospective_precision':'99999 case permutations if rank-association family adopted; store exact integer centered doubled-rank dot products (n10 magnitude<=330) plus all2000casebootstrap values to preserve complete audit within quota. Actual method/input/resource contract must be accepted before true functional tests.','limits':'Functional complement response would support independent target genes, not receptor binding, clinical prediction or formal cell state. Shared upstream and reference-tissue effects remain possible.'}
    (out/'plan_before_reference_complements_or_functional_RNA.json').write_text(json.dumps(plan,indent=2)+'\n')
    wpath=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wpath)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';w=pd.read_csv(wpath,sep='\t');w=w[w.fit_id.eq('K2_seed20261018_source_adjusted')&w['view'].isin(['R_FOXP3_capture','E_effector_capture'])&w.weight.ne(0)];exclude=set(w.gene.str.upper())
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');read=lambda key:np.asarray(read_elem(z[key])).astype(str);source_order=read('uns/source_order').tolist();measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,source_order.index('GSE176021')];symbol=read('var/gene_symbol');canonical=read('var/canonical_gene_identity');identity=read('var/feature_identity_status');alias=read('var/source_gene_aliases');mapping=defaultdict(set)
    for i in np.flatnonzero(measured):
        if identity[i]=='ambiguous_source_feature' or canonical[i].startswith('AMBIGUOUS::'):continue
        names={symbol[i].upper(),canonical[i].upper()}
        if identity[i]=='stable_Ensembl_ID':names.update(g.upper() for g in alias[i].split('|'))
        for name in names:mapping[name].add(int(i))
    details=[];summary=[]
    for r in records:
        p=refs/(r['name']+'.json');assert sha(p)==r['sha256'];full=sorted(set(g.upper() for g in json.loads(p.read_text())[r['name']]['geneSymbols']));removed=sorted(set(full)&exclude);remaining=sorted(set(full)-exclude);mapped={g:next(iter(mapping[g])) for g in remaining if len(mapping[g])==1};retained=len(remaining)/len(full);coverage=len(mapped)/len(remaining) if remaining else 0.;eligible=retained>=.8 and coverage>=.8 and len(mapped)>=20;reason='REFERENCE_METADATA_ELIGIBLE' if eligible else 'PROGRAMME_EXCLUSION_RETAINS_BELOW80PCT' if retained<.8 else 'SOURCE_MEASURED_BELOW80PCT' if coverage<.8 else 'FEWER_THAN20_MEASURED_COMPLEMENT_GENES';summary.append({'reference':r['name'],'full_reference_gene_count':len(full),'excluded_predictor_gene_count':len(removed),'remaining_gene_count':len(remaining),'retained_fraction':retained,'unique_measured_complement_genes':len(mapped),'source_measured_fraction':coverage,'eligible':eligible,'status':reason});details.append({'reference':r['name'],'source_reference_sha256':r['sha256'],'full_reference_genes':full,'removed_frozen_predictor_genes':removed,'remaining_reference_genes':remaining,'unique_source_mapping':mapped,'unmeasured_or_ambiguous_genes':sorted(set(remaining)-set(mapped))})
    pd.DataFrame(summary).to_csv(out/'private_all50_functional_reference_qualification.tsv',sep='\t',index=False);(out/'private_all50_full_and_disjoint_reference_gene_definitions.json').write_text(json.dumps(details,indent=2)+'\n');np.savez_compressed(out/'private_all_source_measured_rank_universe.npz',native_columns=np.flatnonzero(measured).astype(np.int32));files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    d={'status':'COMPLETE_ALL50_DISJOINT_FUNCTION_REFERENCE_QUALIFICATION_NO_RNA_SCORES_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'references':50,'metadata_eligible_references':sum(x['eligible'] for x in summary),'potential_slots':400,'current_pyucell_version':version,'new_functional_scores_or_association':False,'new_gene_fit_or_download':False}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json'),'metadata_eligible_references':d['metadata_eligible_references'],'pyucell_version':version}))


if __name__=='__main__':main()
