"""Independent pandas grouped means for every author RNA-composition sidecar unit."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_author_RNA_composition_20261016_v1'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_author_RNA_composition_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_author_RNA_composition.json').read_text())
    assert plan['code_sha256']==receipt['source_sha256']
    for f,h in terminal['files'].items():assert sha(out/f)==h
    native=r/'NSCLC_full_native_gene_units_20261016_v1';source=r/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz';assert sha(source)==plan['input_sha256']['metadata']
    units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t').sort_values('unit_index');idx=np.load(native/'private_cell_unit_index.npy',allow_pickle=False)
    assert sha(native/'private_cell_unit_index.npy')==plan['input_sha256']['unit_index']
    md=pd.read_csv(source,usecols=['sampleID','sub_cell_type','total_counts','total_counts_mt','pct_counts_mt','total_counts_rb','pct_counts_rb']);assert len(idx)==len(md)
    assert (md.sampleID.to_numpy()==units.sampleID.to_numpy()[idx]).all() and (md.sub_cell_type.to_numpy()==units.sub_cell_type.to_numpy()[idx]).all()
    got=pd.read_csv(out/'private_all_native_unit_RNA_composition.tsv',sep='\t');assert np.array_equal(got.unit_index,np.arange(len(units)))
    grouped=md.assign(unit_index=idx).groupby('unit_index',sort=True);assert np.array_equal(grouped.size(),got.n_cells) and np.array_equal(got.n_cells,units.n_cells)
    maximum=0.
    for short,field,component in [('mito','pct_counts_mt','total_counts_mt'),('ribo','pct_counts_rb','total_counts_rb')]:
        v=pd.to_numeric(md[field],errors='coerce').to_numpy(float);finite=np.isfinite(v);count=pd.Series(finite.astype(int)).groupby(idx,sort=True).sum().to_numpy();complete=count==got.n_cells.to_numpy()
        expected=pd.Series(np.where(finite,v,np.nan)).groupby(idx,sort=True).mean().to_numpy();expected[~complete]=np.nan
        observed=got[f'mean_author_pct_{short}'].to_numpy(float);assert np.allclose(observed,expected,rtol=1e-12,atol=1e-10,equal_nan=True)
        maximum=max(maximum,float(np.nanmax(np.abs(observed-expected))));assert np.array_equal(got[f'finite_{short}_cells'],count) and np.array_equal(got[f'complete_{short}'],complete)
        raw=pd.to_numeric(md[component],errors='coerce').to_numpy(float);total=md.total_counts.to_numpy(float);paired=finite&np.isfinite(raw)
        assert np.max(np.abs(v[paired]-raw[paired]/total[paired]*100))<=plan['percent_formula_tolerance']
    proof=dict(status='PASS_ALL_SOURCE_CELL_ORDER_UNIT_KEYS_COUNTS_AUTHOR_PERCENT_FORMULA_AND_ALL_NATIVE_UNIT_GROUPED_MEANS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),units_checked=len(got),maximum_mean_difference=maximum,new_association_P=False)
    p=out/'private_independent_QC_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
