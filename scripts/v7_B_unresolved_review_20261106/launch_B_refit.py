import os,sys,json,hashlib,subprocess,psutil
from datetime import datetime,timezone
r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1';name='B_unresolved'
assert json.loads((r/'stage_c_QC_views/independent_input_audit.json').read_text())['status']=='PASS'
assert not (r/('stage_d_'+name+'_refit')).exists()
rs=Path(sys.executable).parents[2]/'scrna_r/bin/Rscript';body=r/'scripts/B_lineage_refit_v1.R'
p=subprocess.run([str(rs),'--vanilla','-e','invisible(parse(file=commandArgs(TRUE)[1]))',str(body)],capture_output=True,text=True);assert p.returncode==0,p.stderr
sup='''import os,subprocess,json
from pathlib import Path
from datetime import datetime,timezone
r=Path(os.environ['B_REVIEW_RUN']);name='B_unresolved';start=datetime.now(timezone.utc).isoformat()
p=subprocess.run([os.environ['B_RSCRIPT'],str(r/'scripts/B_lineage_refit_v1.R'),str(r/'stage_c_QC_views'/name/'refit_config.json')],stdin=subprocess.DEVNULL)
(r/'launch'/(name+'_refit_terminal.json')).open('x').write(json.dumps({'UTC_start':start,'UTC_end':datetime.now(timezone.utc).isoformat(),'exit_code':p.returncode,'lineage':name,'complete_exists':(r/('stage_d_'+name+'_refit/complete.json')).exists()},indent=2))
'''
sp=r/'launch/supervise_B_refit.py';sp.open('x').write(sup)
env=os.environ.copy();env.update(PROJECT_ROOT=str(PROJECT_ROOT),B_REVIEW_RUN=str(r),B_RSCRIPT=str(rs),OPENBLAS_NUM_THREADS='16',OMP_NUM_THREADS='16',MKL_NUM_THREADS='16')
with (r/'launch/B_unresolved_refit_stdout_stderr.log').open('xb') as log:
 p=subprocess.Popen([sys.executable,'-B',str(sp)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
out={'UTC':datetime.now(timezone.utc).isoformat(),'lineage':name,'supervisor_pid':p.pid,'create_time':psutil.Process(p.pid).create_time(),'source_R_SHA256':hashlib.sha256(body.read_bytes()).hexdigest(),'threads':16,'input_audit':'PASS','HVG_IG_variable_globin_exclusions':True,'full_RNA_retained':True};(r/'launch/B_unresolved_refit_launch.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
