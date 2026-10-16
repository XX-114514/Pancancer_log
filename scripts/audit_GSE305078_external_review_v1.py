"""Audit all private display rows and original resample geometry, no new P."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from PIL import Image


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);src=r/'GSE305078_external_ligand_programme_core_results_20261016_v1';out=r/'GSE305078_external_core_private_review_20261016_v1'
    report=json.loads((out/'private_terminal.json').read_text());assert report['input_audit_sha256']==sha(src/'private_independent_audit.json')
    for n,h in report['files'].items():assert sha(out/n)==h
    t=pd.read_csv(src/'private_all676_external_core_results.tsv',sep='\t');cfg=json.loads((r/'external_ligand_programme_core_followup_config_20261016_v1.json').read_text())
    forest=pd.read_csv(out/'private_all_estimated_forest_data.tsv',sep='\t');expected=t[t.wild_P.notna()].sort_values(['ligand','readout','model_id']).reset_index(drop=True)
    pd.testing.assert_frame_equal(forest,expected,check_dtype=False)
    c=pd.read_csv(out/'private_all676_coverage.tsv',sep='\t');counts=t.value_counts(['stage','extension_family','readout','status']).sort_index()
    actual=c.set_index(['stage','extension_family','readout','status']).planned_slots.sort_index();assert actual.index.equals(counts.index) and np.array_equal(actual,counts)
    geom=pd.read_csv(out/'private_all_estimated_bootstrap_geometry.tsv',sep='\t').set_index('model_id');checked=0
    for i,row in enumerate(t.itertuples()):
        if not np.isfinite(row.wild_P):continue
        p=src/'private_resamples'/(row.model_id+'.npz');b=np.load(p,allow_pickle=False)['bootstrap_rho'];rng=np.random.default_rng(cfg['seed']+i*10+1)
        counts=np.array([np.count_nonzero(np.bincount(rng.integers(row.library_units,size=row.library_units),minlength=row.library_units)) for _ in range(2000)])
        z=geom.loc[row.model_id];v=b[np.isfinite(b)]
        values=[counts.min(),np.median(counts),np.mean(counts<=row.effective_rank+2),np.mean(np.abs(v)>=1-1e-10),row.bootstrap_high-row.bootstrap_low]
        actual=[z.unique_units_minimum,z.unique_units_median,z.bootstrap_fraction_unique_le_main_rank_plus2,z.finite_bootstrap_fraction_at_abs_rho1,z.original_CI_width]
        assert np.allclose(values,actual,atol=1e-12) and np.isfinite(b).sum()==z.bootstrap_valid and z.no_draws_removed
        checked+=1
    png=out/'private_all_estimable_forest.png'
    with Image.open(png) as im:im.verify()
    with Image.open(png) as im:dimensions=list(im.size)
    assert dimensions==[1500,1350] and checked==len(forest)==len(geom)
    proof=dict(status='PASS_ALL676_COVERAGE_ALL_ESTIMATED_FOREST_ROWS_AND_BOOTSTRAP_GEOMETRY_IMAGE_DECODABLE_VISUAL_REVIEW_SEPARATE',finished_UTC=datetime.now(timezone.utc).isoformat(),
               audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),estimated_rows_checked=checked,forest_png_sha256=sha(png),image_dimensions=dimensions,new_inference=False)
    p=out/'private_independent_display_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
