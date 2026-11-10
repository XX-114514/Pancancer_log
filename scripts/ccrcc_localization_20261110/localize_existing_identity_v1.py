"""Decompose existing ccRCC spatial evidence; never recalculate identities.

All regional tables and coordinates remain private. Public tables aggregate
across blocks; images rasterize cells rather than exporting coordinate rows.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys, json, hashlib, math
import numpy as np
import pandas as pd
from scipy.special import gammaln

SEND = ['Blood_EC', 'Mural', 'Lymphatic']
RECV = ['B', 'T', 'NK']
AXES = [('CCL19', 'CCR7'), ('CXCL12', 'CXCR4'), ('CCL21', 'CCR7')]
GENES = ['CCL19', 'CXCL12', 'CCR7', 'CXCR4']
IDENTITY_SHA = 'd65cbcb6f8508f2598934d9de3c95a072fde4a53ba2871a13689728a8e262bde'


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def strata(xy, depth, tile, identities=None):
    keys = np.column_stack([np.floor(xy / tile).astype(int), np.floor(np.log2(depth)).astype(int)])
    if identities is not None:
        keys = np.column_stack([keys, identities])
    _, inverse = np.unique(keys, axis=0, return_inverse=True)
    return inverse


def missing_probability(n, successes, draws):
    """P(no selected success) in a uniform sample without replacement."""
    n, successes, draws = np.broadcast_arrays(n, successes, draws)
    out = np.zeros(n.shape, dtype=float)
    z = (draws <= n - successes)
    out[z] = np.exp(gammaln(n[z]-successes[z]+1)-gammaln(n[z]-successes[z]-draws[z]+1)
                    +gammaln(n[z]-draws[z]+1)-gammaln(n[z]+1))
    return np.clip(out, 0, 1)


def shared_expectation(owner, receiver_strata, cell_codes, all_strata, n_ec):
    """Exact expectation of EC having >=1 B AND >=1 NK under original null."""
    ng = int(all_strata.max()) + 1
    size = np.bincount(all_strata, minlength=ng)
    nb = np.bincount(all_strata[cell_codes == 0], minlength=ng)
    nn = np.bincount(all_strata[cell_codes == 2], minlength=ng)
    if not len(owner): return np.zeros(n_ec)
    key, k = np.unique(np.column_stack([owner, receiver_strata]), axis=0, return_counts=True)
    ec, st = key.T
    products = []
    for hits in [nb, nn, nb + nn]:
        factors = missing_probability(size[st], hits[st], k)
        product = np.ones(n_ec)
        np.multiply.at(product, ec, factors)
        products.append(product)
    return np.clip(1-products[0]-products[1]+products[2], 0, 1)


def region_summary(ob, ex, den):
    """Fixed-null contribution deletion, not refitted-null inference or patient CI."""
    good = den > 0
    stable = den >= 20
    delta = ob - ex
    total_den = den.sum()
    answer = {'observed': float(ob.sum()), 'exact_expected': float(ex.sum()),
              'difference': float(delta.sum()),
              'difference_pp': float(100*delta.sum()/total_den) if total_den else None,
              'blocks_with_opportunities': int(good.sum()),
              'blocks_at_least20_opportunities': int(stable.sum()),
              'positive_blocks_at_least20': int((stable & (delta > 0)).sum()),
              'positive_fraction_at_least20': float((delta[stable] > 0).mean()) if stable.any() else None,
              'max_positive_block_share_of_all_positive_excess': float(max(0, delta.max())/np.maximum(delta,0).sum()) if np.maximum(delta,0).sum() else None}
    remaining = total_den-den[good]
    effects = np.divide(100*(delta.sum()-delta[good]),remaining,out=np.full(remaining.shape,np.nan),where=remaining>0)
    answer.update({'fixed_null_leave_one_block_out_min_pp': float(np.nanmin(effects)) if np.isfinite(effects).any() else None,
                   'fixed_null_leave_one_block_out_max_pp': float(np.nanmax(effects)) if np.isfinite(effects).any() else None,
                   'positive_after_every_single_block_deletion': bool(np.all(effects > 0)) if len(effects) else None})
    return answer


def main(root, out):
    out.mkdir(exist_ok=False)
    prior = root/'current_identity_forward_v1'
    identity_file = root/'spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz'
    assert digest(identity_file) == IDENTITY_SHA
    d = pd.read_csv(identity_file, sep='\t', index_col=0)
    xy = d[['x_centroid','y_centroid']].to_numpy()
    label = d.identity_candidate.astype(str).to_numpy()
    qc = d.QC_keep.to_numpy(bool)
    depth = d.rna_total.to_numpy()
    assert qc.sum() == 358841 and np.isfinite(xy).all()
    saved = np.load(prior/'communication_inputs_PRIVATE.npz')
    inc = saved['included']; values = saved['expression']
    assert len(inc) == 151833
    meta = pd.read_csv(prior/'imaging/ccRCC_current_identity/metadata_PRIVATE.tsv',sep='\t',index_col=0)
    assert np.array_equal(d.index[inc].astype(str),meta.index.astype(str))
    assert np.array_equal(label[inc],meta.broad.to_numpy())
    es, er, dist, pair = [saved[k] for k in ['es','er','distance','pair_code']]
    assert np.all(label[inc[es]] != label[inc[er]])
    block_xy, blocks = np.unique(np.floor(xy/500).astype(int),axis=0,return_inverse=True)
    nblock = len(block_xy)
    active = np.bincount(blocks[qc],minlength=nblock)
    rows, region_rows = [], []
    prior_lr = pd.read_csv(prior/'spatial_expression_communication_complete108.tsv',sep='\t')
    for tile in [250,500]:
        ids = pd.Categorical(label[inc],categories=SEND+RECV).codes
        st = strata(xy[inc],depth[inc],tile,ids)
        n = np.bincount(st)
        means = np.column_stack([np.bincount(st,weights=values[:,j])/n for j in range(4)])
        for rad in [50,100]:
            for ai,(ligand,receptor) in enumerate(AXES):
                for si,sender in enumerate(SEND):
                    for ri,receiver in enumerate(RECV):
                        prefix = dict(sender=sender,receiver=receiver,ligand=ligand,receptor=receptor,radius_um=rad,tile_um=tile)
                        q = (dist<=rad)&(pair==si*3+ri)
                        br = blocks[inc[er[q]]]
                        den = np.bincount(br,minlength=nblock)
                        if ai==2:
                            rows.append({**prefix,'status':'MISSING_PANEL_COMPONENT','spatial_edges':int(q.sum()),'observed':None,'exact_expected':None})
                            continue
                        li,rr = GENES.index(ligand),GENES.index(receptor)
                        positive = values[es[q],li] & values[er[q],rr]
                        expected = means[st[es[q]],li]*means[st[er[q]],rr]
                        ob = np.bincount(br,weights=positive,minlength=nblock)
                        ex = np.bincount(br,weights=expected,minlength=nblock)
                        old = prior_lr[(prior_lr.sender==sender)&(prior_lr.receiver==receiver)&(prior_lr.ligand==ligand)&(prior_lr.radius_um==rad)&(prior_lr.tile_um==tile)].iloc[0]
                        assert int(ob.sum())==int(old.observed_dualpositive_edges) and int(den.sum())==old.spatial_edges
                        rows.append({**prefix,'status':'MEASURED_EXACT_EXPECTATION_DESCRIPTIVE','spatial_edges':int(q.sum()),**region_summary(ob,ex,den),'previous_permutation_null_mean':float(old.null_mean_dualpositive_edges),'previous_BH_q':float(old.BH_q_complete108),'previous_BY_q':float(old.BY_q_complete108)})
                        for b in np.flatnonzero(den):
                            region_rows.append({**prefix,'block_x':int(block_xy[b,0]),'block_y':int(block_xy[b,1]),'opportunities':int(den[b]),'observed':float(ob[b]),'expected':float(ex[b]),'excess':float(ob[b]-ex[b]),'QC_cells_in_block':int(active[b])})
    lr = pd.DataFrame(rows); assert len(lr)==108
    lr.to_csv(out/'LR108_exact_expectation_block_robustness.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(region_rows).to_csv(out/'LR_regional_PRIVATE.tsv.gz',sep='\t',index=False)
    graph = np.load(prior/'graph_PRIVATE.npz')
    imm = graph['immune_indices']; codes = pd.Categorical(label[imm],categories=RECV).codes
    ec = np.flatnonzero(qc & (label=='Blood_EC'))
    ec_map = np.full(len(d),-1);ec_map[ec] = np.arange(len(ec))
    eg = ec_map[graph['sender']]>=0
    owner = ec_map[graph['sender'][eg]]
    ir = graph['receiver_row'][eg];dd = graph['distance'][eg]
    comp_rows,comp_regional=[],[]
    oldcomp = pd.read_csv(prior/'CM03_existing_identity_spatial12.tsv',sep='\t')
    for tile in [250,500]:
        st = strata(xy[imm],depth[imm],tile)
        sizes = np.bincount(st)
        for rad in [50,100]:
            near = dd<=rad
            for relation,code in [('EC_B_edges',0),('EC_NK_edges',2),('same_EC_B_and_NK',None)]:
                key=dict(radius_um=rad,tile_um=tile,relation=relation)
                if code is not None:
                    p = np.bincount(st,weights=codes==code)/sizes
                    br = blocks[imm[ir[near]]]
                    ob = np.bincount(br,weights=codes[ir[near]]==code,minlength=nblock)
                    ex = np.bincount(br,weights=p[st[ir[near]]],minlength=nblock)
                    den = np.bincount(br,minlength=nblock)
                    opportunity = 'all_EC_to_immune_edges'
                else:
                    nb=np.bincount(owner[near & (codes[ir]==0)],minlength=len(ec))
                    nn=np.bincount(owner[near & (codes[ir]==2)],minlength=len(ec))
                    expected=shared_expectation(owner[near],st[ir[near]],codes,st,len(ec))
                    br=blocks[ec]
                    ob=np.bincount(br,weights=(nb>0)&(nn>0),minlength=nblock)
                    ex=np.bincount(br,weights=expected,minlength=nblock)
                    den=np.bincount(br,minlength=nblock)
                    opportunity='all_Blood_EC_cells'
                old=oldcomp[(oldcomp.radius_um==rad)&(oldcomp.tile_um==tile)&(oldcomp.relation==relation)].iloc[0]
                assert int(ob.sum())==old.observed_count
                comp_rows.append({**key,'opportunity_unit':opportunity,**region_summary(ob,ex,den),'previous_permutation_null_mean':float(old.null_mean),'previous_BH_q':float(old.BH_q_complete12),'previous_BY_q':float(old.BY_q_complete12)})
                for b in np.flatnonzero(den):
                    comp_regional.append({**key,'block_x':int(block_xy[b,0]),'block_y':int(block_xy[b,1]),'opportunities':int(den[b]),'observed':float(ob[b]),'expected':float(ex[b]),'excess':float(ob[b]-ex[b])})
    comp=pd.DataFrame(comp_rows);assert len(comp)==12
    comp.to_csv(out/'CM12_exact_expectation_block_robustness.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(comp_regional).to_csv(out/'CM_regional_PRIVATE.tsv.gz',sep='\t',index=False)
    # Fine-grained grid data remain private. Public figures contain raster summaries only.
    raster_step=100
    origin=np.floor(xy.min(axis=0)/raster_step)*raster_step
    grid=np.floor((xy-origin)/raster_step).astype(int)
    shape=grid.max(axis=0)+1;flat=grid[:,1]*shape[0]+grid[:,0]
    def raster(indices,weights=None):
        return np.bincount(flat[indices],weights=weights,minlength=int(np.prod(shape))).reshape(shape[::-1])
    arrays={k:raster(np.flatnonzero(qc & (label==k))) for k in SEND+RECV}
    arrays['QC']=raster(np.flatnonzero(qc))
    for ident in SEND:
        for gene in GENES[:2]:
            z=label[inc]==ident
            arrays[ident+'_'+gene]=raster(inc[z],values[z,GENES.index(gene)].astype(int))
    for ident in RECV:
        for gene in GENES[2:]:
            z=label[inc]==ident
            arrays[ident+'_'+gene]=raster(inc[z],values[z,GENES.index(gene)].astype(int))
    np.savez_compressed(out/'raster_PRIVATE.npz',origin=origin,step=raster_step,**arrays)
    assert digest(identity_file)==IDENTITY_SHA
    receipt={'status':'COMPUTED_EXISTING_IDENTITY_LOCALIZATION_AUDIT_PENDING','UTC':datetime.now(timezone.utc).isoformat(),'identity_sha256':IDENTITY_SHA,'original_QC_cells':int(qc.sum()),'six_original_identity_cells':len(inc),'composition_slots':len(comp),'communication_slots':len(lr),'new_P_values':0,'identity_changed':False,'coordinates_rescaled':False,'shuffle_model_changed':False,'expectation':'exact marginal means for distinct-identity LR; multivariate hypergeometric inclusion-exclusion for shared EC','localization_block_um':500,'raster_um':100,'pathology_image_in_downloaded_members':False,'no_patient_replication_claim':True}
    (out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
