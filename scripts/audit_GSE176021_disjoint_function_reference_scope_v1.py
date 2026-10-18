"""All50 complement gene definitions/native identities independently reconstructed."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);source=root/'GSE176021_disjoint_function_reference_scope_20261019_v1';out=root/'GSE176021_disjoint_function_reference_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);t=json.loads((source/'private_terminal.json').read_text())
    for n,v in t['files'].items():assert sha(source/n)==v['sha256']
    refs=Path(__file__).parent/'hallmark50_public_reference_v3';records=json.loads((refs/'reference_provenance.json').read_text())['records'];w=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');block={str(r.gene).upper() for r in w.itertuples(index=False) if r.fit_id=='K2_seed20261018_source_adjusted' and r.view in ['R_FOXP3_capture','E_effector_capture'] and r.weight!=0};z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');read=lambda key:np.asarray(read_elem(z[key])).astype(str);order=read('uns/source_order').tolist();measure=np.asarray(read_elem(z['varm/measured_in_source']),bool)[:,order.index('GSE176021')];sym=read('var/gene_symbol');canon=read('var/canonical_gene_identity');typ=read('var/feature_identity_status');alias=read('var/source_gene_aliases');definitions=json.loads((source/'private_all50_full_and_disjoint_reference_gene_definitions.json').read_text());summary=pd.read_csv(source/'private_all50_functional_reference_qualification.tsv',sep='\t').set_index('reference');assert len(definitions)==len(summary)==len(records)==50
    universe=np.load(source/'private_all_source_measured_rank_universe.npz');assert np.array_equal(universe['native_columns'],np.flatnonzero(measure))
    requested=set(g for d in definitions for g in d['full_reference_genes']);hits={g:set() for g in requested}
    for i in np.flatnonzero(measure):
        if typ[i]=='ambiguous_source_feature' or canon[i].startswith('AMBIGUOUS::'):continue
        ns={sym[i].upper(),canon[i].upper()}
        if typ[i]=='stable_Ensembl_ID':ns|={x.upper() for x in alias[i].split('|')}
        for g in ns&requested:hits[g].add(int(i))
    for ref in records:
        p=refs/(ref['name']+'.json');assert sha(p)==ref['sha256'];full={str(g).upper() for g in json.loads(p.read_text())[ref['name']]['geneSymbols']};remain=full-block;removed=full&block;measured={g:next(iter(hits[g])) for g in remain if len(hits[g])==1};d=next(x for x in definitions if x['reference']==ref['name']);assert set(d['full_reference_genes'])==full and set(d['removed_frozen_predictor_genes'])==removed and set(d['remaining_reference_genes'])==remain and d['unique_source_mapping']==measured and set(d['unmeasured_or_ambiguous_genes'])==remain-set(measured);assert not set(measured)&block;r=summary.loc[ref['name']];assert r.full_reference_gene_count==len(full) and r.excluded_predictor_gene_count==len(removed) and r.remaining_gene_count==len(remain) and r.unique_measured_complement_genes==len(measured);assert abs(r.retained_fraction-len(remain)/len(full))<1e-12 and abs(r.source_measured_fraction-len(measured)/len(remain))<1e-12;eligible=len(remain)/len(full)>=.8 and len(measured)/len(remain)>=.8 and len(measured)>=20;assert r.eligible==eligible
    d={'status':'ACCEPTED_ALL50_COMPLEMENT_DEFINITIONS_AND_SOURCE_IDENTITIES_NO_FUNCTIONAL_RNA','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(source/'private_terminal.json'),'references':50,'metadata_eligible_references':int(summary.eligible.sum()),'functional_scores_or_P':False,'runtime_scoring_qualification_pending':True};p=out/'private_terminal.json';p.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(p),'metadata_eligible_references':d['metadata_eligible_references']}))


if __name__=='__main__':main()
