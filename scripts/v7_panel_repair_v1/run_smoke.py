"""Run a synthetic full-matrix workflow and independent arithmetic audit."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rscript',required=True);ap.add_argument('--registry',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);out=a.out.resolve();code=Path(__file__).resolve().parent
    shutil.copyfile(a.registry,out/'registry.json');start=datetime.now(timezone.utc).isoformat();stages=[]
    def run(stage,args):
        with (out/(stage+'.log')).open('wb') as log:p=subprocess.run([a.rscript,*map(str,args)],stdout=log,stderr=subprocess.STDOUT)
        stages.append({'stage':stage,'exit_code':p.returncode})
        (out/'terminal.json').write_text(json.dumps({'start_UTC':start,'updated_UTC':datetime.now(timezone.utc).isoformat(),'stages':stages,'complete':stage=='verify' and p.returncode==0},indent=2)+'\n')
        if p.returncode:raise RuntimeError(stage+' failed; see private log')
    run('fixture',[code/'smoke_fixture.R',out])
    run('workflow',[code/'lineage_workflow.R',out/'input.rds',out/'config.json',out/'result'])
    run('verify',[code/'verify_smoke.R',out])
    print((out/'independent_arithmetic_audit.json').read_text())
if __name__=='__main__':main()
