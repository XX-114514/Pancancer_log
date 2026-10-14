"""Independent metadata selection replay, with no expression or outcome use."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,math
import numpy as np,pandas as pd
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
I=R/'all_CM_state_native_input_preflight_20261014_v1';D=R/'all_lineage_state_discovery_row_selection_20261014_v1'
plan=json.loads(Path(os.environ['PANCANCER_STATE_DISCOVERY_ROW_PLAN']).read_text());iv=json.loads((D/'inventory.json').read_text())
assert sha(I/'inventory.json')==plan['input_inventory_sha256']
u=pd.read_csv(I/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False)
f=np.load(I/'private_canonical_row_selection.npz');z=np.load(D/'private_discovery_row_indices.npz')
idx=z['selected_input_indices'];assert len(idx)==100766 and (np.diff(idx)>0).all()
assert np.array_equal(z['canonical_global_rows'],f['global_rows'][idx]) and np.array_equal(z['unit_indices'],f['unit_index'][idx])
m=pd.read_csv(D/'private_discovery_metadata.tsv.gz',sep='\t',keep_default_na=False)
expected=u.iloc[f['unit_index'][idx]].reset_index(drop=True)
for c in u.columns:assert np.array_equal(m[c].to_numpy(),expected[c].to_numpy()),c
assert np.array_equal(m.core_position,z['canonical_global_rows'])
s=pd.read_csv(D/'all18_major_discovery_scope.tsv',sep='\t');q=pd.read_csv(D/'all_selected_source_major_discovery_scope.tsv',sep='\t')
genes=json.loads((D/'fixed_common_training_genes_by_major.json').read_text());replayed=[];checks=[]
allmaj=sorted(u.major.unique());assert len(allmaj)==18
for mi,major in enumerate(allmaj):
    v=u[u.major.eq(major)];records=[]
    for sid,group in v.groupby('gse_id',sort=True):
        records.append((sid,int(group.n_cells.sum()),group.sample_uid.nunique(),group.effective_resampling_block.nunique()))
    eligible=sorted([r for r in records if r[1]>=plan['min_source_cells'] and r[3]>=plan['min_source_blocks']],key=lambda r:(-r[3],-r[2],-r[1],r[0]))[:plan['max_sources_per_major']]
    sids=sorted(r[0] for r in eligible);actual=q[q.major.eq(major)];assert sorted(actual.gse_id)==sids
    cap=math.ceil(plan['max_cells_per_major']/max(len(sids),1));gsets=[];n=0
    for sj,sid in enumerate(sids):
        allowed_units=np.flatnonzero(u.major.eq(major)&u.gse_id.eq(sid));which=np.flatnonzero(np.isin(f['unit_index'],allowed_units))
        blocks=u.effective_resampling_block.to_numpy()[f['unit_index'][which]];rng=np.random.default_rng(plan['seed']+mi*100+sj)
        queues=[rng.permutation(which[blocks==b]) for b in sorted(set(blocks))]
        # Interleave the independently constructed queues by depth and block order.
        take=[];depth=0;target=min(len(which),cap)
        while len(take)<target:
            take.extend(int(a[depth]) for a in queues if len(a)>depth)
            depth+=1
        take=take[:target];replayed.extend(take);n+=len(take)
        observed=np.flatnonzero(m.major.eq(major)&m.gse_id.eq(sid));assert len(observed)==target
        r=actual[actual.gse_id.eq(sid)].iloc[0];assert r.available_cells==len(which) and r.training_cells==target
        selected_blocks=m.iloc[observed].effective_resampling_block.value_counts()
        capacities=pd.Series(blocks).value_counts();rounds=0
        while rounds<int(capacities.max()) and sum(min(c,rounds+1) for c in capacities)<=target:rounds+=1
        for b,c in capacities.items():assert min(c,rounds)<=selected_blocks.get(b,0)<=min(c,rounds+1)
        gp=AR/'20261010_uniform_expression_annotation_v1/by_gse'/sid/'genes.tsv';gsets.append(set(pd.read_csv(gp,sep='\t').gene))
    common=sorted(set.intersection(*gsets));assert common==genes[major]
    row=s[s.major.eq(major)].iloc[0];assert row.selected_cells==n and row.common_training_genes==len(common)
    assert row.source_limited==(len(sids)<2) and row.model_training_eligible==(n>=plan['min_training_cells'] and len(common)>=plan['min_common_genes'])
    assert row.available_cells==int(v.n_cells.sum()) and row.available_sources==len(records)
    checks.append(dict(major=major,selected_cells=n,selected_sources=len(sids),common_genes=len(common)))
assert np.array_equal(np.sort(replayed),idx)
for name,h in iv['public_output_sha256'].items():assert sha(D/name)==h
assert sha(D/'private_discovery_row_indices.npz')==iv['private_row_sha256']
out=dict(status='PASS_ALL100766_DISCOVERY_ROWS_ALL18_MAJORS_RANK_CAP_BLOCK_REPLAY_AND_GENE_AXES',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_inventory_sha256=sha(D/'inventory.json'),audit_code_sha256=sha(__file__),selected_cells=len(idx),all_labels_and_canonical_rows_exact=True,metadata_only_fixed_selection=True,patient_independence_not_assumed=True,source_limited_majors=3,checks=checks,formal_states=False)
O=R/'all_lineage_state_discovery_row_audit_20261014_v1';assert not O.exists();O.mkdir();(O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
