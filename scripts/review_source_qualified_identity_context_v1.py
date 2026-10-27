"""Join source qualification to native group identity and carried QC evidence.

Extra contextual genes are unused by the original manual marker panels, but
may occur in reference models or integration features. They are NOT independent
validation. This review changes no identity, QC, membership or clinical status.
"""
import argparse, collections, datetime, gc, hashlib, json, resource, sys, time
from pathlib import Path
import numpy as np

CONTEXT_PANELS = {
    'T_receptor_context': ['TRBC1', 'TRBC2', 'LCK', 'LAT', 'ZAP70'],
    'NK_receptor_context': ['NCR1', 'KLRF1', 'KLRD1', 'GNLY'],
    'B_receptor_context': ['CD19', 'CD22', 'CD79B', 'PAX5'],
    'plasma_secretory_context': ['XBP1', 'SSR4', 'FKBP11', 'SEC11C'],
    'myeloid_context': ['TYROBP', 'SPI1', 'CD68', 'MS4A7'],
    'endothelial_context': ['CDH5', 'CLDN5', 'ESAM', 'ERG'],
    'lymphatic_context': ['LYVE1', 'MMRN1', 'CCL21', 'PROX1'],
    'fibroblast_context': ['PDGFRA', 'COL3A1', 'COL6A1', 'COL6A2'],
    'contractile_context': ['MYH11', 'CNN1', 'TAGLN', 'ACTA2'],
    'epithelial_context': ['KRT7', 'KRT14', 'KRT5', 'SFN'],
    'glial_context': ['MOG', 'MOBP', 'MAG', 'CNP'],
}

def load(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, d): Path(p).write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
def hist(v):
    names, n = np.unique(np.asarray(v).astype(str), return_counts=True)
    return dict(zip(names.tolist(), map(int,n)))

def extra_context(n, cp, detected, genes, excluded_genes):
    """Retain per-gene quantities and missingness, never a new type classifier."""
    genes=list(genes); n=np.asarray(n); index={g:i for i,g in enumerate(genes)}
    if len(index)!=len(genes) or n.ndim!=1 or (n<=0).any():
        raise ValueError('Unique measured genes and positive exact-group sizes required')
    if cp.shape!=(len(n),len(genes)) or detected.shape!=cp.shape:
        raise ValueError('Native group/gene axes differ')
    panels={k:[g for g in v if g not in set(excluded_genes)] for k,v in CONTEXT_PANELS.items()}
    selected=list(dict.fromkeys(g for v in panels.values() for g in v))
    values={}
    for g in selected:
        if g not in index: values[g]=None; continue
        c=cp[:,index[g]]; d=detected[:,index[g]]
        if not np.isfinite(c).all() or not np.isfinite(d).all() or (c<0).any() or (d<0).any() or (d>n).any():
            raise ValueError('Invalid native sufficient statistic')
        values[g]={'mean_CP10K':c/n,'detection_fraction':d/n}
    coverage={k:{'requested':v,'measured':[g for g in v if g in index],
                 'not_measured':[g for g in v if g not in index]} for k,v in panels.items()}
    return coverage,values

def source_status(atlas_source, core, manifest):
    """Apply current and inherited per-cell overrides, not parent white lists."""
    ds=atlas_source['source']; rows=atlas_source['rows']; ua=atlas_source['ua']
    latest=manifest['new_source_overrides']; older=load(manifest['base_v4_manifest']['path'])
    # Original mixed author samples must stay split at their native cell axis.
    original=older['bindings']['author_bridge_registry']; assert sha(original['path'])==original['sha256']
    inherited_source=load(original['path'])['source']
    all_overrides={inherited_source:{'path':older['override_arrays_path'],
                                    'sha256':older['override_arrays_sha256']}}
    all_overrides.update(latest)
    if ds in all_overrides:
        r=all_overrides[ds];p=Path(r['path']);assert sha(p)==r['sha256']
        with np.load(p,allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_row'],rows)
            # v4 and v5 work-unit dictionary positions differ; status is the
            # explicit inherited scientific decision, not the old position code.
            status=f['source_qualification_code'].copy()
    else:
        by={u['qualification_unit_key_sha256']:u for u in core}
        decode={'排除':0,'来源资格通过':1,'待核实':2}
        code=np.array([decode[by[u['native_unit_key_sha256']]['source_qualification']]
                       for u in atlas_source['um']['unit_records']],np.int8)
        status=code[ua['native_expression_unit_code']]
    if status.shape!=(len(rows),) or not np.isin(status,[0,1,2]).all():
        raise ValueError('Invalid source status axis')
    return status

def main(root):
    started=time.monotonic();cfg=load(root/'config/spec_v1.json')
    for b in cfg['bindings'].values(): assert sha(b['path'])==b['sha256']
    joined=Path(cfg['joined']);sys.path.insert(0,str(joined/'scripts'))
    from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
    from native_annotation_review_v1 import REVIEW_GENES
    from audit_full_source_local_graphs_v1 import index_at
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    atlas=JoinedNativeAtlas(joined/'outputs/joined_native_atlas_delivery_v1.json')
    marker=Path(cfg['markers']);mi=index_at(marker,'full_source_local_marker_evidence_index_v1.json')
    view=load(cfg['source_view']);core=load(view['work_units']['path'])['units']
    assert sha(view['work_units']['path'])==view['work_units']['sha256']
    before=load(view['base_v4_manifest']['path']);assert sha(view['base_v4_manifest']['path'])==view['base_v4_manifest']['sha256']
    totals=collections.Counter();summary=[];index=[]
    qfields=['preintegration_cell_eligible','preintegration_release_eligible','downstream_release_eligible',
             'residual_doublet_flag','annotation_exclude_from_primary','prior_QC_status','QC_prior',
             'doublet_status','predicted_doublet','scrublet_success','strict_major_analysis_eligible',
             'authoritative_l1_downstream_eligible']
    for ds in atlas.sources:
        storage_gate(root,cfg)
        s=atlas.source(ds);ca=s['ca'];rows=s['rows'];fit=ca['fit_pool_candidate'].astype(bool)
        status=source_status(s,core,view);qualified=status==1;active=qualified&fit
        cm=s['cm']['candidate_dictionaries'];maj=np.asarray(cm['consensus_major'],object)[ca['consensus_major_code']]
        fine=np.asarray(cm['consensus_celltype'],object)[ca['consensus_celltype_code']]
        reason=np.asarray(cm['coarse_evidence'],object)[ca['coarse_evidence_code']]
        freq=np.asarray(cm['fine_evidence'],object)[ca['fine_evidence_code']]
        available=set(atlas.native.field_names(ds));qvalues={f:atlas.native.field(ds,f) for f in qfields if f in available}
        mpath=marker/'outputs'/mi[ds]['metadata'];assert sha(mpath)==mi[ds]['metadata_sha256'];m=load(mpath)
        stat=marker/'outputs'/m['group_native_statistics'];assert sha(stat)==m['group_native_statistics_sha256']
        with np.load(stat,allow_pickle=False) as f:
            n=f['n'];cp=f['cp_sum'];detected=f['detected']
        coverage,context=extra_context(n,cp,detected,m['native_measured_gene_names'],REVIEW_GENES)
        group_code=ca['guarded_source_subset_code'];assert (group_code[fit]>=0).all()
        np.testing.assert_array_equal(np.bincount(group_code[fit],minlength=len(n)),n)
        qn=np.bincount(group_code[active],minlength=len(n))
        pairs=np.unique(np.column_stack([group_code[fit],ca['global_subcluster_number'][fit]]),axis=0)
        np.testing.assert_array_equal(pairs[:,0],np.arange(len(n)))
        global_by_code=pairs[:,1]
        records=[]
        for j,descriptor in enumerate(m['guarded_subsets']):
            assert descriptor['guarded_subset_code']==j and descriptor['n_subset_fit_rows']==int(n[j])
            g=atlas.groups[int(global_by_code[j])]
            records.append({'source':ds,'guarded_subset_code':j,'global_subcluster_number':int(global_by_code[j]),
                'stable_group_sha256':g['stable_group_sha256'],'native_fit_group_rows':int(n[j]),
                'source_qualified_fit_rows':int(qn[j]),'current_candidate_celltype':descriptor['guarded_celltype_candidate'],
                'subcluster_display':g['subcluster'],'own_marker':g['marker'],
                'marker_evidence_type':g['marker_evidence_type'],'group_expression_evidence_status':g['group_expression_evidence_status'],
                'n_comparable_expression_units_not_patients':g['n_comparable_expression_units_not_patients'],
                'own_top_markers':[x['gene'] for x in m['groups'][j]['own_group_markers']],
                'additional_native_gene_evidence':{gene:None if v is None else {k:float(a[j]) for k,a in v.items()} for gene,v in context.items()},
                'same_cell_coexpression_for_additional_genes':'not_available_from_group_sums',
                'additional_genes_may_be_in_models_not_independent_validation':True})
        # Preserve every unresolved reason; no denominator is formed by silently
        # dropping unresolved or non-fitted archive cells.
        r={'source':ds,'n_archive_rows':len(rows),'source_status_counts':hist(status),
           'n_source_qualified_rows':int(qualified.sum()),'n_source_qualified_fit_candidates_not_final_QC':int(active.sum()),
           'source_qualified_coarse_candidates':hist(maj[qualified]),'source_qualified_fine_candidates':hist(fine[qualified]),
           'source_qualified_coarse_evidence':hist(reason[qualified]),'source_qualified_fine_evidence':hist(freq[qualified]),
           'source_qualified_fit_coarse_candidates':hist(maj[active]),
           'carried_QC_flag_distributions':{f:hist(v[qualified]) for f,v in qvalues.items()},
           'carried_QC_fields_absent':[f for f in qfields if f not in available],
           'n_groups':len(records),'n_groups_with_source_qualified_fit_rows':int((qn>0).sum()),
           'n_qualified_rows_without_fitted_group':int((qualified&~fit).sum()),
           'carried_flag_names_do_not_certify_final_QC':True}
        out=root/'outputs'/(ds+'_source_qualified_identity_context_v1.json')
        save(out,{'summary':r,'context_panel_coverage':coverage,'groups':records,'input_native_statistics_SHA256':sha(stat),
                  'manual_panel_genes_excluded_from_context':True,'no_relabeling_or_final_QC_call':True})
        index.append({'source':ds,'path':str(out),'SHA256':sha(out)})
        summary.append(r);totals.update(n_sources=1,n_archive_rows=len(rows),n_source_qualified_rows=int(qualified.sum()),
            n_source_qualified_fit_candidates_not_final_QC=int(active.sum()),n_groups=len(records),
            n_groups_with_source_qualified_fit_rows=int((qn>0).sum()))
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'n_completed_sources':len(summary),'source':ds})
        del s,ca,cp,detected,n,context,records,m,qvalues;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
    expected_archive=sum(u['n_archive_rows'] for u in core)
    expected_qualified=sum(u['n_archive_rows'] for u in core if u['source_qualification']=='来源资格通过')
    assert totals['n_archive_rows']==expected_archive and totals['n_source_qualified_rows']==expected_qualified and totals['n_groups']==len(atlas.groups)
    assert sum(sum(r['source_qualified_coarse_candidates'].values()) for r in summary)==expected_qualified
    save(root/'outputs/source_qualified_identity_context_index_v1.json',index)
    save(root/'outputs/source_qualified_identity_context_summary_v1.json',{'totals':dict(totals),'sources':summary,
        'original_manual_marker_genes_removed_from_context':REVIEW_GENES,'clinical_outcomes_not_used':True,
        'identity_QC_and_specimen_gates_not_upgraded':True,'source_qualification_is_not_final_CM_eligibility':True})
    save(root/'provenance/terminal_v1.json',{'status':'SOURCE_QUALIFIED_ALL_ARCHIVE_IDENTITY_CONTEXT_AND_QC_EVIDENCE_REVIEW_COMPLETE',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),
        'summary_SHA256':sha(root/'outputs/source_qualified_identity_context_summary_v1.json'),
        'index_SHA256':sha(root/'outputs/source_qualified_identity_context_index_v1.json'),
        'native_count_matrix_read':False,'new_annotation_or_QC_labels':False,'wall_seconds':time.monotonic()-started,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try: main(r)
    except BaseException as e:
        import traceback
        save(r/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
