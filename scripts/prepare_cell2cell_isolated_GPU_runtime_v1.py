"""Small isolated official wheels, reusing an existing compatible CUDA runtime."""
import json
import os
import sys
from datetime import datetime,timezone
from pathlib import Path


PACKAGES={'cell2cell':'0.8.0','tensorly':'0.9.0','gseapy':'1.0.3','kneed':'0.8.5','statannotations':'0.7.2','xlrd':'2.0.2'}


def worker():
    import hashlib
    import importlib.metadata
    import inspect
    import io
    import stat
    import urllib.request
    import zipfile
    from packaging.tags import sys_tags
    from packaging.utils import parse_wheel_filename
    import numpy as np
    import torch
    out=Path(os.environ['PANCANCER_CELL2CELL_RUNTIME']);site=out/'site_packages';assert site.exists() and not list(site.iterdir())
    assert np.__version__=='1.26.4' and torch.cuda.is_available()
    assert float((torch.arange(4,device='cuda')**2).sum().item())==14.
    compatible={tag:i for i,tag in enumerate(sys_tags())};wheels=[];metadata=[]
    for name,version in PACKAGES.items():
        url='https://pypi.org/pypi/'+name+'/'+version+'/json'
        with urllib.request.urlopen(url,timeout=30) as response:raw=response.read(2*1024**2+1)
        assert len(raw)<=2*1024**2;a=json.loads(raw);assert a['info']['version']==version
        choices=[]
        for record in a['urls']:
            if record['packagetype']!='bdist_wheel':continue
            _,_,_,tags=parse_wheel_filename(record['filename']);match=[compatible[t] for t in tags if t in compatible]
            if match:choices.append((min(match),record['size'],record['filename'],record))
        assert choices,'NO_COMPATIBLE_OFFICIAL_BINARY_WHEEL'
        record=sorted(choices,key=lambda x:x[:3])[0][-1];assert record['url'].startswith('https://files.pythonhosted.org/') and record['size']<=16*1024**2
        with urllib.request.urlopen(record['url'],timeout=40) as response:data=response.read(16*1024**2+1)
        assert len(data)==record['size'] and hashlib.sha256(data).hexdigest()==record['digests']['sha256']
        wheel=zipfile.ZipFile(io.BytesIO(data));unpacked=0
        for info in wheel.infolist():
            p=Path(info.filename);assert not p.is_absolute() and '..' not in p.parts and not stat.S_ISLNK(info.external_attr>>16)
            assert not any(part.endswith('.data') for part in p.parts),'WHEEL_DATA_LAYOUT_REQUIRES_STANDARD_INSTALLER'
            unpacked+=info.file_size
        wheels.append((wheel,data));metadata.append({'name':name,'version':version,'official_metadata_url':url,'wheel_url':record['url'],'wheel_name':record['filename'],'wheel_sha256':record['digests']['sha256'],'compressed_bytes':len(data),'uncompressed_bytes':unpacked,'requirements':a['info']['requires_dist']})
    expected=sum(d['uncompressed_bytes'] for d in metadata);assert expected+4*1024**2<48*1024**2
    (out/'private_official_binary_dependency_manifest.json').write_text(json.dumps({'observed_UTC':datetime.now(timezone.utc).isoformat(),'packages':metadata,'predicted_extracted_bytes':expected,'no_wheels_or_CUDA_download_persisted':True},indent=2)+'\n')
    for wheel,_ in wheels:wheel.extractall(site)
    sys.path.insert(0,str(site))
    import cell2cell as c2c
    import tensorly as tl
    import gseapy
    tl.set_backend('pytorch')
    names=list(PACKAGES)+['numpy','scipy','pandas','torch','scanpy','umap-learn','statsmodels','seaborn']
    versions={name:importlib.metadata.version(name) for name in names}
    for name,version in PACKAGES.items():assert versions[name]==version
    from tensorly.decomposition import non_negative_parafac
    receipt={'status':'AVAILABLE_OFFICIAL_CELL2CELL_GPU_RUNTIME_SYNTHETIC_CP_GATE_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'python':sys.version.split()[0],'versions':versions,'cuda':torch.version.cuda,'GPU_tiny_fixture_sum':14,'tensor_backend':tl.get_backend(),'APIs':{'PrebuiltTensor':str(inspect.signature(c2c.tensor.PrebuiltTensor)),'run_tensor_cell2cell_pipeline':str(inspect.signature(c2c.analysis.run_tensor_cell2cell_pipeline)),'non_negative_parafac':str(inspect.signature(non_negative_parafac))},'real_expression_or_model_used':False,'isolated_site_package_files':sum(p.is_file() for p in site.rglob('*'))}
    (out/'private_worker_runtime_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    import subprocess
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,48*1024**2);out=root/'Tensor_cell2cell_GPU_runtime_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);(out/'site_packages').mkdir()
    python=Path(sys.executable).resolve().parents[2]/'scvi/bin/python';assert python.exists()
    (out/'plan_before_isolated_official_wheels.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'packages':PACKAGES,'existing_runtime':'existingPython3.11CUDA12.6Torch;NumPy1.26.4compatibility required before install','scope':'Official compatible binary wheels only; SHA and safezip/uncompressed size inspected before extraction. Isolated newsite directory,no basepackage changes,no CUDA/torch/numpy downloads. Original wheels RAMonly; no deletion/overwrite of existing artifacts.','gate':'Software import/API and tinyGPU arithmetic only. Full missingmask/CP/heldoutprojection synthetic test required before real communication data or clinical inference.'},indent=2)+'\n')
    env=os.environ.copy();env.update(PANCANCER_CELL2CELL_RUNTIME=str(out),PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(out/'private_numba_cache'),MPLCONFIGDIR=str(out/'private_matplotlib_cache'),XDG_CACHE_HOME=str(out/'private_external_cache'),MPLBACKEND='Agg')
    with (out/'private_install_worker_stdout.log').open('xb') as stdout,(out/'private_install_worker_stderr.log').open('xb') as stderr:r=subprocess.run([str(python),'-u',str(Path(__file__).resolve()),'--worker'],env=env,stdout=stdout,stderr=stderr,timeout=420)
    assert r.returncode==0,'ISOLATED_BINARY_RUNTIME_WORKER_FAILED'
    d=json.loads((out/'private_worker_runtime_receipt.json').read_text());files={str(p.relative_to(out)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};actual=sum(p.stat().st_size for p in out.rglob('*') if p.is_file());assert actual<48*1024**2
    receipt={'status':'AVAILABLE_OFFICIAL_CELL2CELL_GPU_RUNTIME_FULL_SYNTHETIC_GATE_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'actual_all_runtime_bytes':actual,'existing_CUDA_software_reused':True,'new_real_expression_or_model':False,'software_versions':d['versions']}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json'),'bytes':actual}))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:worker()
    else:main()
