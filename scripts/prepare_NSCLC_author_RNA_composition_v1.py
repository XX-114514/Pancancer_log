"""Preserve author mitochondrial/ribosomal RNA percentages as native-unit sidecars."""
import hashlib,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

FIELDS=['total_counts','total_counts_mt','pct_counts_mt','total_counts_rb','pct_counts_rb']


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);native=r/'NSCLC_full_native_gene_units_20261016_v1'
    paths=dict(metadata=r/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',unit_index=native/'private_cell_unit_index.npy',units=native/'private_all_native_units_and_quality.tsv',native_audit=r/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json',clinical_scope=r/'NSCLC_all_case_clinical_scope_20261016_v1/private_terminal.json')
    old=json.loads((native/'private_terminal.json').read_text());audit=json.loads(paths['native_audit'].read_text())
    assert audit['status'].startswith('PASS_ALL_OUTPUTS_OLD_OVERLAP') and sha(native/'private_terminal.json')==audit['producer_terminal_sha256']
    hashes={k:sha(p) for k,p in paths.items()}
    assert hashes['metadata']==json.loads(paths['clinical_scope'].read_text())['source_metadata_sha256']
    for k,name in [('unit_index','private_cell_unit_index.npy'),('units','private_all_native_units_and_quality.tsv')]:assert hashes[k]==old['outputs'][name]['sha256']
    out=r/'NSCLC_author_RNA_composition_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    df=subprocess.run(['df','-P',str(r.resolve())],capture_output=True,text=True,timeout=20);quota=subprocess.run(['quota','-g','-v','-w','USER002'],capture_output=True,text=True,timeout=20)
    assert df.returncode==quota.returncode==0;fs=df.stdout.splitlines()[-1].split()[0];qr=[x.split() for x in quota.stdout.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1
    remaining=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert remaining>500_000_000_000+1_073_741_824+268_435_456
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256=hashes,
              scope='Every original source cell and every native unit, no clinical or effect filtering.',
              measurement='Arithmetic mean of original author per-cell pct_counts_mt and pct_counts_rb; not pooled-count fraction. Verify count/total percent formula within 0.001 percentage point.',
              missing='Any nonfinite author percentage in a unit makes its mean unavailable; preserve valid zero. Do not exclude cases to obtain association estimates.',
              interpretation='RNA composition sensitivity, may represent technical quality or real metabolic/ribosomal state. Not automatic artifact removal or a new state label.',
              future_scope='All212 original programme/subtype slots for composition associations and all10600 functional slots with both extra covariates; no new P/BH or old family modification.',
              percent_formula_tolerance=1e-3,quota_free_bytes=remaining,output_cap_bytes=268_435_456)
    (out/'plan_before_author_RNA_composition.json').write_text(json.dumps(plan,indent=2)+'\n')
    md=pd.read_csv(paths['metadata'],usecols=['sampleID','sub_cell_type']+FIELDS);u=np.load(paths['unit_index'],allow_pickle=False);units=pd.read_csv(paths['units'],sep='\t').sort_values('unit_index')
    assert len(md)==len(u)==old['native_cells'] and np.array_equal(units.unit_index,np.arange(len(units)))
    assert u.min()==0 and u.max()==len(units)-1
    assert np.array_equal(md.sampleID.to_numpy(),units.sampleID.to_numpy()[u]) and np.array_equal(md.sub_cell_type.to_numpy(),units.sub_cell_type.to_numpy()[u])
    counts=np.bincount(u,minlength=len(units));assert np.array_equal(counts,units.n_cells.to_numpy(int))
    total=pd.to_numeric(md.total_counts,errors='coerce').to_numpy(float);assert np.isfinite(total).all() and (total>0).all()
    result=pd.DataFrame(dict(unit_index=np.arange(len(units)),n_cells=counts));formula_errors={};missing={}
    for short,component,percent in [('mito','total_counts_mt','pct_counts_mt'),('ribo','total_counts_rb','pct_counts_rb')]:
        v=pd.to_numeric(md[percent],errors='coerce').to_numpy(float);component_v=pd.to_numeric(md[component],errors='coerce').to_numpy(float);valid=np.isfinite(v)
        assert ((v[valid]>=0)&(v[valid]<=100)).all();paired=valid&np.isfinite(component_v)
        assert (component_v[paired]>=0).all() and (component_v[paired]<=total[paired]).all()
        error=float(np.max(np.abs(v[paired]-100*component_v[paired]/total[paired]))) if paired.any() else np.nan
        assert np.isfinite(error) and error<=plan['percent_formula_tolerance'];formula_errors[short]=error
        valid_counts=np.bincount(u,weights=valid.astype(int),minlength=len(units)).astype(int);sums=np.bincount(u,weights=np.where(valid,v,0),minlength=len(units));means=sums/counts
        complete=valid_counts==counts;means[~complete]=np.nan
        result[f'mean_author_pct_{short}']=means;result[f'finite_{short}_cells']=valid_counts;result[f'complete_{short}']=complete;missing[short]=int((~complete).sum())
    result.to_csv(out/'private_all_native_unit_RNA_composition.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_AUTHOR_RNA_COMPOSITION_ALL_NATIVE_UNITS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),native_cells=len(md),native_units=len(units),maximum_percent_formula_errors=formula_errors,incomplete_units=missing,
                  files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_association_P=False,new_expression_copy=False)
    assert sum(p.stat().st_size for p in out.iterdir())<plan['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
