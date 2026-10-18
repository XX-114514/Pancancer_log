"""Freeze all three lineage models, all50 complements and bounded native assay."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


STRESS={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'}
VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    assert sha(root/'NSCLC_original_author_lineage_audit_20261019_v1/private_terminal.json')=='ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d'
    assert sha(root/'NSCLC_immune_RNA_MCP_native_axis_audit_20261019_v1/private_terminal.json')=='575bb2a84376ddd6557882ca47f0f8045d2be7c5d08854082482b0da099759e4'
    assert sha(root/'immune_RNA_MCP_all_atlas_native_measurement_audit_20261019_v1/private_terminal.json')=='f7a6db048cfccac2a98de84023ce5048d4faa595353d2589767718d0e661cfc0'
    gene_path=root/'GSE243013_genes.csv.gz';assert sha(gene_path)=='ac325c9ce6bfaf2270eb931e7d947bf31d1fb151b14cb166ed9f6c8f4b34a104'
    genes=pd.read_csv(gene_path,keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str);assert len(genes)==len(set(genes))==31831
    source=root/'NSCLC_original_author_lineage_qualification_20261019_v1';pool=pd.read_csv(source/'private_all729_author_lineage_pool_positions.tsv',sep='\t');member=dict(np.load(source/'private_all_author_lineage_native_memberships.npz'))
    selected=member['native_row'][member['selected']];assert len(selected)==68300 and len(set(selected))==len(selected)
    out=root/'NSCLC_author_lineage_assay_scope_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    native=root/'NSCLC_full_native_gene_units_20261016_v1';detected=np.load(native/'private_detected_genes.npy',mmap_mode='r');selected_nnz=int(np.sum(detected[selected],dtype=np.int64));memory_needed=selected_nnz*80+6*1024**3
    bound=len(selected)*(50*20+4*8+4*8+4+2+8+4)+729*50*4*8+8*1024**2;assert bound<96*1024**2
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'author_metadata_acceptance_sha256':'ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d','scope':'All691author-eligible pools/729registered positions,all619247author lineage cells/68300fixedhashselected. Both MCPs, all3source author lineages and full/stress. Myeloid remains broad17author categories, not macrophage-only; CD8union includes authorMAIT. No formal state/global annotation mutation.','primary_measurement':'Sameoriginalhash2026101802 up100 perpool; full native logCP10K and fixed weights/scales/coefficients. Primary nuisance fraction/logN retains original marker/major capture count and allimmune denominator; ownauthor-population adjustment separately labelled. Unknown additive source constant remains omitted. Fullauthor-cell programme means also measured.','functional_reference':'All50 original public Hallmark sets. Exclude union of active genes in ALL3views/BOTHMCPs and any physical native overlap, not the old two-view complement. >=80% original retention,>=80% unique native measurement,>=20mappedgenes; failed references remain missing slots. Common newthree-view complement across all roles.','functional_scores':'Existing pyucell0.8.0 CPU/average ties/skip,full31831nativefeatures rankbackground;1500primary,1000/2000andmeanlogCP10Ksensitivities. Only same68300hashselected cells. No raw-cell count matrix written to disk; full selected sparse matrix stays RAM until scorecompletion.','native_audit':'Existing6.6GBcompressedmatrix localonly,no download. Stream full2,010,550,708sortedinteger coordinates,verify compressedSHA/gzipCRC and allcellRNA/detected. Fullindependent native replay/allnew selectedprogram/quality/function values and fullpopulationmeans required before effect review. Every oldoverlap primaryscore checked,Myeloid anchor unchanged.','prospective_output_cap_bytes':96*1024**2,'conservative_output_bound_bytes':bound,'selected_native_nnz_for_RAM_guard':selected_nnz,'RAM_guard_required_bytes':memory_needed,'inference':'No true new clinical/functionassociation here. Complete newsource families and mathematical/resource/design gates must be separately frozen. Potential900functional slots are metadata scope only, not registered/executed tests. Oldfinite6NSCLC/400GSE176021families remainclosed.','record_units':'Authorrecords only; no reviewed independent patient aliases. Twoidentical nonoutcomeprofiles not merged/declared samepatient. Source postneoadjuvant protocol not perrecord newtiming annotation.'}
    (out/'plan_before_all_author_native_RNA.json').write_text(json.dumps(plan,indent=2)+'\n')
    wp=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wp)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';w=pd.read_csv(wp,sep='\t');w=w[w.fit_id.eq('K2_seed20261018_source_adjusted')&w['view'].isin(VIEWS)]
    params=pd.read_csv(root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2/private_all12_frozen_models.tsv',sep='\t').set_index(['view','program','variant']);old=pd.read_csv(root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1/private_original_discovery_nuisance_coefficients.tsv',sep='\t').set_index(['view','program','term']).coefficient
    gm=pd.read_csv(root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1/private_all_active_gene_native_axis.tsv',sep='\t').set_index(['view','program','gene']).native_column
    beta=np.zeros((3,len(genes),4));centers=np.zeros((3,4));coef=np.zeros((3,4,6));mu=np.zeros((3,4));sd=np.zeros((3,4));parameters=[];active=set()
    train=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    for vi,view in enumerate(VIEWS):
        z=np.load(train/(view+'_frozen_projection.npz'));mu[vi]=z['quality_mean'];sd[vi]=z['quality_sd']
        for pi,pr in enumerate(['MCP1','MCP2']):
            ww=w[w['view'].eq(view)&w.program.eq(pr)]
            for si,var in enumerate(['full','fixed_stress_complement']):
                use=ww if si==0 else ww[~ww.gene.str.upper().isin(STRESS)];columns=np.array([gm.loc[(view,pr,g)] for g in use.gene],int);assert (columns>=0).all() and len(set(columns))==len(columns)
                b=use.weight.to_numpy()/use.discovery_sd.to_numpy();k=2*pi+si;beta[vi,columns,k]=b;centers[vi,k]=float(b@use.discovery_mean.to_numpy());r=params.loc[(view,pr,var)];terms=['intercept','cellQ','geneQ','sourceCC','fraction','logN'];coef[vi,k]=r[terms].to_numpy(float);assert abs(centers[vi,k]-r.center)<1e-8 and len(use)==r.active_genes
                if si==0:
                    for ti,term in enumerate(terms):assert abs(coef[vi,k,ti]-old.loc[(view,pr,term)])<1e-8
                    active.update(columns[use.weight.ne(0)].tolist())
                parameters.append({'view':view,'program':pr,'variant':var,'active_gene_positions':len(use),'complete_native_axis':True})
    assert (sd>0).all();np.savez_compressed(out/'private_all12_frozen_native_models.npz',native_beta=beta,center=centers,coef=coef,quality_mean=mu,quality_sd=sd,views=np.array(VIEWS))
    pd.DataFrame(parameters).to_csv(out/'private_all12_model_native_qualification.tsv',sep='\t',index=False)
    references=Path(__file__).parent/'hallmark50_public_reference_v3';manifest=references/'reference_provenance.json';prior=json.loads((root/'immune_RNA_MCP_functional_description_20261018_v1/plan_before_complete_functional_description.json').read_text());assert sha(manifest)==prior['public_reference_manifest_sha256'];records=json.loads(manifest.read_text())['records'];exclude=set(w[w.weight.ne(0)].gene.str.upper());lookup={g:i for i,g in enumerate(genes)};details=[];qual=[]
    for record in records:
        p=references/(record['name']+'.json');assert sha(p)==record['sha256'];full=sorted(set(g.upper() for g in json.loads(p.read_text())[record['name']]['geneSymbols']));removed=sorted(g for g in full if g in exclude or lookup.get(g,-1) in active);remain=sorted(set(full)-set(removed));mapped={g:lookup[g] for g in remain if g in lookup};assert not set(mapped.values())&active
        retain=len(remain)/len(full);cover=len(mapped)/len(remain) if remain else 0.;eligible=retain>=.8 and cover>=.8 and len(mapped)>=20;status='ELIGIBLE_ALL3VIEW_PREDICTOR_DISJOINT_COMPLEMENT' if eligible else 'RETAINED_REFERENCE_BELOW80PCT' if retain<.8 else 'SOURCE_MEASUREMENT_BELOW80PCT' if cover<.8 else 'FEWER_THAN20_MEASURED_GENES'
        details.append({'reference':record['name'],'original_public_reference_sha256':record['sha256'],'full_reference_genes':full,'excluded_all3view_frozen_genes':removed,'remaining_genes':remain,'unique_source_mapping':mapped,'eligible':eligible,'status':status});qual.append({'reference':record['name'],'original_genes':len(full),'excluded_predictor_genes':len(removed),'remaining_genes':len(remain),'unique_source_genes':len(mapped),'retained_fraction':retain,'source_measured_fraction':cover,'eligible':eligible,'status':status})
    assert len(details)==50
    (out/'private_all50_three_view_complement_definitions.json').write_text(json.dumps(details,indent=2)+'\n');pd.DataFrame(qual).to_csv(out/'private_all50_reference_qualification.tsv',sep='\t',index=False)
    (out/'private_original_input_SHA_manifest.json').write_text(json.dumps({'metadata_source_sha256':sha(root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'),'source_genes_sha256':sha(gene_path),'weights_sha256':sha(wp),'public_reference_manifest_sha256':sha(manifest),'author_metadata_terminal_sha256':sha(source/'private_terminal.json'),'download_receipt_sha256':sha(root/'count_matrix_download_v1/download_state.json'),'old_official_UCell_worker_sha256':sha(Path(__file__).parent/'measure_GSE176021_disjoint_function_RNA_v1.py')},indent=2)+'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt={'status':'COMPLETE_ALL3_AUTHOR_LINEAGE_NATIVE_ASSAY_MODELS50_REFERENCES_AND_STORAGE_SCOPE_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'selected_cells':len(selected),'author_eligible_pools':int(pool.eligible.sum()),'source_features':len(genes),'model_variants':12,'all_reference_positions':50,'metadata_eligible_references':sum(x['eligible'] for x in qual),'new_expression_or_effects':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
