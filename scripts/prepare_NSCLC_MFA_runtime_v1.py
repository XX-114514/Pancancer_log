"""Install one verified official wheel privately; synthetic data only in smoke test."""
import hashlib,importlib.metadata,json,os,subprocess,sys,urllib.request
from datetime import datetime,timezone
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_MFA_runtime_20261016_v1';assert not out.exists();fs=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);rows=[s.split() for s in quota.splitlines() if s.split() and s.split()[0]==fs];assert len(rows)==1;free=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456;out.mkdir(mode=0o700)
    versions={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','h5py','anndata','scikit-learn']};metadata=json.load(urllib.request.urlopen('https://pypi.org/pypi/mofapy2/0.7.2/json',timeout=30));wheels=[x for x in metadata['urls'] if x['filename']=='mofapy2-0.7.2-py3-none-any.whl'];assert len(wheels)==1;item=wheels[0];expected='7d1e73e1d189c9ff3bdd529d91201fbdc4dac635b4c9bae6b08a4902563abb49';assert item['digests']['sha256']==expected
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),package='mofapy2',version='0.7.2',wheel_sha256=expected,base_versions=versions,install='Official PyPI verified wheel, private --target library, --no-deps --no-cache-dir; no base environment mutation.',software_download_only=True,new_expression_or_clinical_input_download=False,quota_free_bytes=free)
    (out/'plan_before_install.json').write_text(json.dumps(plan,indent=2)+'\n');wheel=out/item['filename'];wheel.write_bytes(urllib.request.urlopen(item['url'],timeout=30).read());assert sha(wheel)==expected and wheel.stat().st_size==item['size']
    with (out/'private_install.log').open('x') as log:subprocess.run([sys.executable,'-m','pip','install','--no-deps','--no-cache-dir','--no-compile','--target',str(out/'library'),str(wheel)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
    sys.path.insert(0,str(out/'library'));from mofapy2.run.entry_point import entry_point
    rng=np.random.default_rng(2026101731);z=rng.normal(size=(60,3));views=[z@rng.normal(size=(3,120))+.3*rng.normal(size=(60,120)) for _ in range(4)];ent=entry_point();ent.set_data_options(scale_views=False,scale_groups=False,center_groups=True,use_float32=False);ent.set_data_matrix([[x] for x in views],likelihoods=['gaussian']*4,views_names=['synthetic_view'+str(i) for i in range(4)],groups_names=['synthetic_group']);ent.set_model_options(factors=3,spikeslab_factors=False,spikeslab_weights=False,ard_factors=False,ard_weights=True);ent.set_train_options(iter=150,convergence_mode='medium',dropR2=None,seed=2026101731,quiet=True,gpu_mode=False);ent.build();ent.run();expectations=ent.model.getExpectations();assert expectations['Z']['E'].shape==(60,3) and np.isfinite(expectations['Z']['E']).all()
    for w in expectations['W']:assert w['E'].shape==(120,3) and np.isfinite(w['E']).all()
    ent.save(str(out/'synthetic_only_model.hdf5'),save_data=False);assert all(importlib.metadata.version(k)==v for k,v in versions.items());status=dict(status='PASS_PINNED_PRIVATE_MOFAPY2_RUNTIME_AND_SYNTHETIC_FOUR_VIEW_FIT_SAVE',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),wheel_sha256=expected,base_versions_unchanged=True,real_expression_or_clinical_data_used=False,files={p.name:sha(p) for p in out.iterdir() if p.is_file()});(out/'private_terminal.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status))


if __name__=='__main__':main()
