"""Small isolated official wheels, reusing an existing compatible CUDA runtime."""
import json
import os
import sys
from datetime import datetime,timezone
from pathlib import Path


PACKAGES={'requests':'2.32.5','charset-normalizer':'3.4.4','idna':'3.11','urllib3':'2.5.0','certifi':'2025.10.5'}


def worker():
    import hashlib
    import importlib.metadata
    import inspect
    import io
    import stat
    import subprocess
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
            unpacked+=info.file_size
        wheels.append((wheel,data));metadata.append({'name':name,'version':version,'official_metadata_url':url,'wheel_url':record['url'],'wheel_name':record['filename'],'wheel_sha256':record['digests']['sha256'],'compressed_bytes':len(data),'uncompressed_bytes':unpacked,'requirements':a['info']['requires_dist']})
    expected=sum(d['uncompressed_bytes'] for d in metadata);assert expected+sum(d['compressed_bytes'] for d in metadata)+4*1024**2<48*1024**2
    (out/'private_official_binary_dependency_manifest.json').write_text(json.dumps({'observed_UTC':datetime.now(timezone.utc).isoformat(),'packages':metadata,'predicted_extracted_bytes':expected,'only_small_official_wheel_archives_persisted_no_CUDA_download':True},indent=2)+'\n')
    wheel_dir=out/'official_wheels';wheel_dir.mkdir();paths=[]
    for (wheel,data),record in zip(wheels,metadata):
        p=wheel_dir/record['wheel_name'];p.write_bytes(data);paths.append(str(p))
    with (out/'private_standard_installer_stdout.log').open('xb') as so,(out/'private_standard_installer_stderr.log').open('xb') as se:
        result=subprocess.run([sys.executable,'-m','pip','install','--no-index','--no-deps','--no-cache-dir','--target',str(site),*paths],stdout=so,stderr=se,timeout=180)
    assert result.returncode==0,'STANDARD_LOCAL_BINARY_WHEEL_INSTALL_FAILED'
    sys.path.insert(0,str(site))
    import cell2cell as c2c
    import tensorly as tl
    import gseapy
    tl.set_backend('pytorch')
    names=list(PACKAGES)+['cell2cell','tensorly','gseapy','kneed','statannotations','xlrd','numpy','scipy','pandas','torch','scanpy','umap-learn','statsmodels','seaborn']
    versions={name:importlib.metadata.version(name) for name in names}
    for name,version in PACKAGES.items():assert versions[name]==version
    from packaging.requirements import Requirement
    seen=set();todo=['cell2cell','tensorly','gseapy','requests']
    while todo:
        name=todo.pop()
        if name.lower() in seen:continue
        seen.add(name.lower());distribution=importlib.metadata.distribution(name)
        for text in distribution.requires or []:
            requirement=Requirement(text)
            if requirement.marker and not requirement.marker.evaluate({'extra':''}):continue
            version=importlib.metadata.version(requirement.name)
            assert not requirement.specifier or requirement.specifier.contains(version,prereleases=True)
            todo.append(requirement.name)
    from tensorly.decomposition import non_negative_parafac
    receipt={'status':'AVAILABLE_OFFICIAL_CELL2CELL_GPU_RUNTIME_SYNTHETIC_CP_GATE_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'python':sys.version.split()[0],'versions':versions,'cuda':torch.version.cuda,'GPU_tiny_fixture_sum':14,'tensor_backend':tl.get_backend(),'APIs':{'PrebuiltTensor':str(inspect.signature(c2c.tensor.PrebuiltTensor)),'run_tensor_cell2cell_pipeline':str(inspect.signature(c2c.analysis.run_tensor_cell2cell_pipeline)),'non_negative_parafac':str(inspect.signature(non_negative_parafac))},'real_expression_or_model_used':False,'isolated_site_package_files':sum(p.is_file() for p in site.rglob('*'))}
    (out/'private_worker_runtime_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    import subprocess
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,48*1024**2);out=root/'Tensor_cell2cell_GPU_runtime_20261019_v3';assert not out.exists();out.mkdir(mode=0o700);(out/'site_packages').mkdir()
    python=Path(sys.executable).resolve().parents[2]/'scvi/bin/python';assert python.exists()
    (out/'plan_before_isolated_official_wheels.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'packages':PACKAGES,'existing_runtime':'existingPython3.11CUDA12.6Torch;NumPy1.26.4compatibility required before install','scope':'Official compatible binary wheels only; SHA/safezip/total compressed+uncompressed size inspected; standard pip localwheel installer handles legitimate .data layout. Isolated newsite directory,no basepackage changes,no CUDA/torch/numpy downloads. Small SHA-verified original wheels retained; no deletion/overwrite of existing artifacts. V1 legitimate .data layout rejection and V2 successful local-wheel install but missing requests import are preserved. V3 supplements only requests and its four dependencies in a separate layer over unchanged V2. All installed dependency specifications are checked below.','gate':'Software import/API and tinyGPU arithmetic only. Full missingmask/CP/heldoutprojection synthetic test required before real communication data or clinical inference.'},indent=2)+'\n')
    base_site=root/'Tensor_cell2cell_GPU_runtime_20261019_v2/site_packages';assert base_site.exists()
    env=os.environ.copy();env.update(PYTHONPATH=str(base_site),PANCANCER_CELL2CELL_RUNTIME=str(out),PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(out/'private_numba_cache'),MPLCONFIGDIR=str(out/'private_matplotlib_cache'),XDG_CACHE_HOME=str(out/'private_external_cache'),MPLBACKEND='Agg')
    with (out/'private_install_worker_stdout.log').open('xb') as stdout,(out/'private_install_worker_stderr.log').open('xb') as stderr:r=subprocess.run([str(python),'-u',str(Path(__file__).resolve()),'--worker'],env=env,stdout=stdout,stderr=stderr,timeout=420)
    assert r.returncode==0,'ISOLATED_BINARY_RUNTIME_WORKER_FAILED'
    d=json.loads((out/'private_worker_runtime_receipt.json').read_text());files={str(p.relative_to(out)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};actual=sum(p.stat().st_size for p in out.rglob('*') if p.is_file());assert actual<48*1024**2
    receipt={'status':'AVAILABLE_OFFICIAL_CELL2CELL_GPU_RUNTIME_FULL_SYNTHETIC_GATE_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'actual_all_runtime_bytes':actual,'existing_CUDA_software_reused':True,'new_real_expression_or_model':False,'software_versions':d['versions']}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json'),'bytes':actual}))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:worker()
    else:main()
