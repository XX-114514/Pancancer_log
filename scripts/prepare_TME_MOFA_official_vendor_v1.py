"""Install one SHA-pinned official wheel into a task-only vendor directory."""
import argparse,base64,datetime,hashlib,importlib.metadata,json,subprocess,sys
from pathlib import Path

def main(R):
 cfg=json.loads((R/'config/spec_v1.json').read_text());raw=base64.b64decode((R/'config/frozen_official_wheel_b64.txt').read_text());assert hashlib.sha256(raw).hexdigest()==cfg['wheel_SHA256']
 wheel=R/'config'/cfg['wheel_filename']
 with wheel.open('xb') as f:f.write(raw)
 vendor=R/'vendor';assert not vendor.exists()
 command=[sys.executable,'-m','pip','install','--no-deps','--no-index','--disable-pip-version-check','--no-compile','--target',str(vendor),str(wheel)]
 p=subprocess.run(command,capture_output=True,text=True,timeout=180)
 for name,body in [('installer_stdout.txt',p.stdout),('installer_stderr.txt',p.stderr)]:
  with (R/'provenance'/name).open('x') as f:f.write(body)
 assert p.returncode==0,'Task vendor install failed; global environment untouched'
 sys.path.insert(0,str(vendor));import mofapy2
 from mofapy2.run.entry_point import entry_point
 version=importlib.metadata.version('mofapy2');assert version==cfg['version'] and Path(mofapy2.__file__).resolve().is_relative_to(vendor.resolve())
 evidence=dict(schema='official_MOFA_task_vendor_v1',version=version,wheel_filename=wheel.name,wheel_SHA256=cfg['wheel_SHA256'],public_distribution=cfg['public_distribution'],vendor_path=str(vendor),imports_successful=True,no_dependencies_installed=True,global_environment_unmodified=True,scientific_data_not_read=True)
 out=R/'outputs/accepted_official_MOFA_vendor_v1.json'
 with out.open('x') as f:json.dump(evidence,f,indent=2)
 with (R/'provenance/terminal_v1.json').open('x') as f:json.dump(dict(status='OFFICIAL_SHA_PINNED_MOFA_TASK_VENDOR_IMPORT_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(out.read_bytes()).hexdigest()),f,indent=2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  with (R/'provenance/failure_v1.json').open('x') as f:json.dump(dict(traceback=traceback.format_exc()),f,indent=2)
  raise
