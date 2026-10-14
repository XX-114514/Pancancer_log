from pathlib import Path
from datetime import datetime,timezone
import os,sys,json,subprocess,tarfile,hashlib,sysconfig
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'state_model_runtime_20261014_v1';bundle=Path(os.environ['PANCANCER_STATE_SOFTWARE_BUNDLE']);assert not O.exists();O.mkdir()
manifest=R/'state_model_runtime_20261014_v1.execution.json'
def now():return datetime.now(timezone.utc).isoformat()
def write(j):manifest.write_text(json.dumps(j,indent=2)+'\n')
j=dict(status='RUNNING_SOFTWARE_INSTALL',started_UTC=now(),pid=os.getpid(),installer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),software_only=True);write(j);code=1
try:
 resources=O/'resources';resources.mkdir()
 with tarfile.open(bundle) as f:
  for m in f.getmembers():assert not m.name.startswith('/') and '..' not in Path(m.name).parts and not m.issym() and not m.islnk()
  f.extractall(resources)
 receipt=json.loads((resources/'runtime_resource_receipts.json').read_text())
 for item in receipt['source_files']:
  assert hashlib.sha256((resources/'spectra_snapshot'/item['path']).read_bytes()).hexdigest()==item['sha256']
 for item in receipt['pinned_wheels']:assert hashlib.sha256((resources/item['filename']).read_bytes()).hexdigest()==item['sha256']
 base='/data4/01.projects/USER002/01.xiongxin/miniconda3/envs/scvi_py311/bin/python';subprocess.run([base,'-m','venv','--system-site-packages',str(O/'venv')],check=True)
 py=str(O/'venv/bin/python');wheels=[str(resources/v['filename']) for v in receipt['pinned_wheels']]
 subprocess.run([py,'-m','pip','install','--ignore-installed','--no-index','--no-deps',*wheels],check=True)
 source=resources/'spectra_snapshot/src';site=subprocess.check_output([py,'-c','import sysconfig;print(sysconfig.get_paths()["purelib"])'],text=True).strip();(Path(site)/'pancancer_pinned_spectra.pth').write_text(str(source)+'\n')
 test=r"""import torch,numpy as np,pandas as pd,anndata as ad,json,importlib.metadata
 import Spectra,cnmf
 from Spectra.Spectra_gpu import est_spectra
 torch.set_num_threads(2);torch.manual_seed(20261014);np.random.seed(20261014)
 a=ad.AnnData(np.log1p(np.random.poisson(2,(12,8))).astype(np.float32));a.var_names=['G'+str(i) for i in range(8)];a.obs['major']=['T cell']*6+['Macrophage']*6
 d={'global':{'shared':['G0','G1','G2']},'T cell':{'T':['G2','G3','G4']},'Macrophage':{'M':['G4','G5','G6']}}
 model=est_spectra(a,d,L={'global':2,'T cell':2,'Macrophage':2},use_highly_variable=False,cell_type_key='major',filter_sets=False,num_epochs=3,batch_size=6,verbose=False,n_top_vals=3)
 assert np.isfinite(a.obsm['SPECTRA_cell_scores']).all() and torch.cuda.is_available()
 print(json.dumps(dict(Python=__import__('sys').version.split()[0],torch=torch.__version__,CUDA=torch.version.cuda,Spectra_commit='9c079a584ab71ea83e1e0ed471f678d2f2a997a2',cnmf=importlib.metadata.version('cnmf'),pyvis=importlib.metadata.version('pyvis'),synthetic_gpu_Spectra_scores_shape=list(a.obsm['SPECTRA_cell_scores'].shape),no_biological_state_fit=True)))
 """
 test='\n'.join(line[1:] if line.startswith(' ') else line for line in test.splitlines())
 with (O/'smoke.stdout.log').open('x') as output:
  p=subprocess.Popen([py,'-c',test],stdout=output,stderr=subprocess.STDOUT);j['child_pid']=p.pid;write(j);rc=p.wait()
 assert rc==0,'Smoke failed; preserve runtime and logs'
 code=0
except Exception as e:j.update(error_type=type(e).__name__,error=str(e)[:250])
finally:
 j.update(status='TERMINAL_SUCCESS_PINNED_GPU_RUNTIME' if code==0 else 'TERMINAL_FAILED_RUNTIME_PREPARATION',finished_UTC=now(),exit_code=code);write(j)
print(json.dumps(j),flush=True)
raise SystemExit(code)
