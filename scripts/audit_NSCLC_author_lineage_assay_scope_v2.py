"""Independent complete frozen model/reference/source-scope qualification."""
import gzip
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_author_lineage_assay_scope_audit_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'NSCLC_author_lineage_assay_scope_20261019_v1';t=wait(root,out,'prepare_NSCLC_author_lineage_assay_scope_v1',src/'private_terminal.json',{'COMPLETE_ALL3_AUTHOR_LINEAGE_NATIVE_ASSAY_MODELS50_REFERENCES_AND_STORAGE_SCOPE_AUDIT_PENDING'})
    for name,v in t['files'].items():assert sha(src/name)==v['sha256']
    assert sha(root/'NSCLC_original_author_lineage_audit_20261019_v1/private_terminal.json')=='ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d'
    axis=root/'NSCLC_immune_RNA_MCP_native_axis_20261019_v1';assert sha(root/'NSCLC_immune_RNA_MCP_native_axis_audit_20261019_v1/private_terminal.json')=='575bb2a84376ddd6557882ca47f0f8045d2be7c5d08854082482b0da099759e4'
    original_axis=json.loads((axis/'private_terminal.json').read_text())
    for name,v in original_axis['files'].items():assert sha(axis/name)==v['sha256']
    native_map=pd.read_csv(axis/'private_all_active_gene_native_axis.tsv',sep='\t').set_index(['view','program','gene']).native_column
    source_genes=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.strip().str.upper().tolist();lookup={name:i for i,name in enumerate(source_genes)};assert len(lookup)==31831
    w=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');w=w[w.fit_id.eq('K2_seed20261018_source_adjusted')];views=['R_FOXP3_capture','E_effector_capture','Myeloid'];assert set(w['view'])==set(views)
    stress={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'}
    models=dict(np.load(src/'private_all12_frozen_native_models.npz'));assert models['views'].tolist()==views and models['native_beta'].shape==(3,31831,4)
    qualification=pd.read_csv(src/'private_all12_model_native_qualification.tsv',sep='\t').set_index(['view','program','variant']);train=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs';physical=set()
    for vi,view in enumerate(views):
        md=pd.read_csv(train/(view+'_discovery_metadata.tsv'),sep='\t');D=np.column_stack([np.ones(len(md)),md[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]);Q,R=np.linalg.qr(D,mode='reduced');gene_order=Path(str(train/(view+'_discovery'))+'_genes.tsv').read_text().splitlines();idx={g:i for i,g in enumerate(gene_order)}
        with gzip.open(str(train/(view+'_discovery'))+'_tpm.mtx.gz','rb') as f:expression=mmread(f).tocsr()
        frozen=np.load(train/(view+'_frozen_projection.npz'));np.testing.assert_array_equal(models['quality_mean'][vi],frozen['quality_mean']);np.testing.assert_array_equal(models['quality_sd'][vi],frozen['quality_sd'])
        for pi,program in enumerate(['MCP1','MCP2']):
            ww=w[w['view'].eq(view)&w.program.eq(program)]
            for si,variant in enumerate(['full','fixed_stress_complement']):
                use=ww if si==0 else ww[~ww.gene.str.upper().isin(stress)];b=np.zeros(31831);beta=use.weight.to_numpy()/use.discovery_sd.to_numpy();columns=[int(native_map.loc[(view,program,g)]) for g in use.gene];assert len(columns)==len(set(columns));b[columns]=beta
                center=float(np.sum(beta*use.discovery_mean.to_numpy()));y=np.asarray(expression[[idx[g] for g in use.gene]].T.multiply(beta).sum(axis=1)).ravel()-center;coefficient=np.linalg.solve(R,Q.T@y);k=2*pi+si
                np.testing.assert_array_equal(models['native_beta'][vi,:,k],b);np.testing.assert_allclose(models['center'][vi,k],center,atol=1e-10,rtol=0);np.testing.assert_allclose(models['coef'][vi,k],coefficient,atol=1e-8,rtol=0)
                rr=qualification.loc[(view,program,variant)];assert rr.complete_native_axis and rr.active_gene_positions==len(use)
                if si==0:physical.update(c for c,v in zip(columns,use.weight) if v!=0)
    refs=Path(__file__).parent/'hallmark50_public_reference_v3';manifest=json.loads((refs/'reference_provenance.json').read_text());definitions=json.loads((src/'private_all50_three_view_complement_definitions.json').read_text());summary=pd.read_csv(src/'private_all50_reference_qualification.tsv',sep='\t').set_index('reference');full_records={r['name']:r for r in manifest['records']};assert len(definitions)==len(summary)==len(full_records)==50
    excluded=set(w[w.weight.ne(0)].gene.str.upper())
    for d in definitions:
        name=d['reference'];record=full_records[name];p=refs/(name+'.json');assert sha(p)==record['sha256']==d['original_public_reference_sha256'];full=set(g.upper() for g in json.loads(p.read_text())[name]['geneSymbols']);removed={g for g in full if g in excluded or lookup.get(g,-1) in physical};remaining=full-removed;mapped={g:lookup[g] for g in sorted(remaining) if g in lookup};retained=len(remaining)/len(full);coverage=len(mapped)/len(remaining) if remaining else 0.;eligible=retained>=.8 and coverage>=.8 and len(mapped)>=20
        assert d['full_reference_genes']==sorted(full) and d['excluded_all3view_frozen_genes']==sorted(removed) and d['remaining_genes']==sorted(remaining) and d['unique_source_mapping']==mapped and d['eligible']==eligible
        r=summary.loc[name];assert [r.original_genes,r.excluded_predictor_genes,r.remaining_genes,r.unique_source_genes]==[len(full),len(removed),len(remaining),len(mapped)];np.testing.assert_allclose([r.retained_fraction,r.source_measured_fraction],[retained,coverage],rtol=0,atol=1e-12);assert r.eligible==eligible and r.status==d['status'] and not set(mapped.values())&physical
    plan=json.loads((src/'plan_before_all_author_native_RNA.json').read_text());membership=np.load(root/'NSCLC_original_author_lineage_qualification_20261019_v1/private_all_author_lineage_native_memberships.npz');selected=membership['native_row'][membership['selected']];detected=np.load(root/'NSCLC_full_native_gene_units_20261016_v1/private_detected_genes.npy',mmap_mode='r')
    bound=len(selected)*(50*(12+8)+32+32+4+2+8+4)+729*50*4*8+8*1024**2;assert plan['conservative_output_bound_bytes']==bound<plan['prospective_output_cap_bytes']==96*1024**2
    assert plan['selected_native_nnz_for_RAM_guard']==sum(int(v) for v in detected[selected]) and plan['RAM_guard_required_bytes']==plan['selected_native_nnz_for_RAM_guard']*80+6*1024**3
    provenance=json.loads((src/'private_original_input_SHA_manifest.json').read_text());assert provenance['old_official_UCell_worker_sha256']=='a18a83257cb21e3377486e1b720ccf17e371f544dd71de0798a8ec93fc7c8ed0'
    receipt={'status':'ACCEPTED_COMPLETE_ALL12_MODELS_ALL50_THREE_VIEW_COMPLEMENTS_AND_BOUNDED_ASSAY_SCOPE','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'native_model_variants':12,'reference_positions':50,'selected_cells':len(selected),'conservative_native_output_bound_bytes':bound,'no_native_assay_or_new_effects':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
