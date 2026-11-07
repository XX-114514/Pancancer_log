"""Finite two-branch continuation; reuse completed preparation, never overwrite scientific stages."""
import sys,subprocess,json,time,hashlib,concurrent.futures
from pathlib import Path
from datetime import datetime,timezone
r=Path(sys.argv[1]);project=r.parents[1];rb=sys.argv[2]
def now():return datetime.now(timezone.utc).isoformat()
def run(tag,cmd):
 receipt={'stage':tag,'UTC_start':now()};assert not (r/'launch'/(tag+'_terminal.json')).exists()
 with (r/'launch'/(tag+'.log')).open('xb') as f:q=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
 receipt.update(UTC_end=now(),exit_code=q.returncode);(r/'launch'/(tag+'_terminal.json')).open('x').write(json.dumps(receipt));return q.returncode

def py(name,args=()):
 code="import runpy,sys;from pathlib import Path;sys.path.insert(0,"+repr(str(r/'scripts'))+");sys.argv=["+repr(name)+"]+"+repr(list(args))+";runpy.run_path("+repr(str(r/'scripts'/name))+",init_globals={'PROJECT_ROOT':Path("+repr(str(project))+")})";return [sys.executable,'-B','-c',code]
def bulk():
 source=r/'sources/SU2C_exact_V4_panel_v1';ready=False
 for _ in range(240):
  if (source/'exact_panel_derivation_receipt.json').exists() and (source/'Table_S1_Clinical_Annotations.xlsx').exists():
   proof=json.loads((source/'exact_panel_derivation_receipt.json').read_text());f=source/'V4_fixed_gene_TPM_PRIVATE.tsv.gz';assert hashlib.sha256(f.read_bytes()).hexdigest()==proof['selected_TPM_SHA256'];ready=True;break
  time.sleep(15)
 (r/'audit/SU2C_full_input_wait_receipt.json').open('x').write(json.dumps({'UTC':now(),'full_author_derived_exact_panel_SHA_ready':ready,'partial_or34gene_input_used':False}))
 # Even when unavailable, all six requested cohorts remain as explicit unassessable slots.
 status={'fit':run('bulk_analysis_v1',py('run_V4_bulk_v2.py'))}
 if status['fit']==0:
  status['measurement_audit']=run('bulk_measurement_audit_v1',py('audit_V4_bulk_measurement.py'));status['numeric_models_audit']=run('bulk_models_audit_v1',[sys.executable,'-B',str(r/'scripts/audit_V4_bulk_models.py'),'--run',str(r/'bulk/analysis_v1'),'--output',str(r/'audit/bulk_models'),'--expected-slots','2100'])
 return status

def imaging():
 results=[];jobs=json.loads((r/'config/imaging_jobs_PRIVATE.json').read_text())
 for i,j in enumerate(jobs):
  z={'unit_code':j['unit_code'],'dataset':j['dataset'],'prepare':run('imaging_'+j['unit_code']+'_prepare',py('prepare_imaging_V4.py',[str(i)]))}
  if z['prepare']==0:
   z['context']=run('imaging_'+j['unit_code']+'_broad_context',py('measure_imaging_broad_CM_context.py',[str(i)]));z['CellChat']=run('imaging_'+j['unit_code']+'_CellChat',[rb,'--vanilla',str(r/'scripts/run_spatial_CellChat_V4.R'),str(r),j['unit_code']])
   if z['CellChat']==0:z['CellChat_audit']=run('imaging_'+j['unit_code']+'_CellChat_audit',[rb,'--vanilla',str(r/'scripts/audit_spatial_CellChat_V4.R'),str(r),j['unit_code']])
  results.append(z);(r/'launch/imaging_progress_PRIVATE.json').write_text(json.dumps({'UTC':now(),'units':results,'n_requested':len(jobs)},indent=2))
 return results
for stage in ['prepare','signatures','CellChat_DB']:assert json.loads((r/'launch'/(stage+'_terminal.json')).read_text())['exit_code']==0
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 a=pool.submit(bulk);b=pool.submit(imaging);summary={'UTC':now(),'bulk':a.result(),'imaging':b.result()}
(r/'launch/extension_all_terminal_PRIVATE.json').open('x').write(json.dumps(summary,indent=2))
