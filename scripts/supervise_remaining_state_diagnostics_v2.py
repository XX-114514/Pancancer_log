"""Bounded supervisor for fixed jobs, completion records and predeclared audits."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,time,hashlib
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'remaining_state_diagnostic_supervision_20261014_v2';assert not O.exists();O.mkdir();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();sp=json.loads(Path(os.environ['PANCANCER_ALL_LINEAGE_SPECTRA_DIAGNOSTIC_PLAN']).read_text());cn=json.loads(Path(os.environ['PANCANCER_ALL_LINEAGE_CNMF_DIAGNOSTIC_PLAN']).read_text());A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';accepted=json.loads(A.read_text());majors=sorted(r['major'] for r in accepted['files']);assert len(majors)==18;ledger=json.loads((V/'recovery_input_inventory.json').read_text());events=[];done=set();start=datetime.now(timezone.utc).isoformat();env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
def record(kind,scope,status,detail=None):
 e=dict(UTC=datetime.now(timezone.utc).isoformat(),kind=kind,scope=scope,status=status,detail=detail);events.append(e);(O/'progress.json').write_text(json.dumps(dict(status='FIXED_SCOPE_DIAGNOSTIC_SUPERVISION',started_UTC=start,events=events,completed_scopes=len(done)),indent=2)+chr(10));print(json.dumps(e),flush=True)
for role,scope,producer,audit in [('Spectra','lymphoid_T',R/'all_lineage_Spectra_diagnostics_20261014_v1/lymphoid_T/inventory.json',R/'all_lineage_Spectra_diagnostic_audit_20261014_v1/lymphoid_T/independent_audit.json'),('cNMF','B cell',R/'all_lineage_cNMF_diagnostics_20261014_v1/major_00/inventory.json',R/'all_lineage_cNMF_diagnostic_audit_20261014_v2/major_00/independent_audit.json')]:
 a=json.loads(audit.read_text());assert a['status'].startswith('PASS_') and a['producer_inventory_sha256']==sha(producer);prefix='lymphoid_T/' if role=='Spectra' else 'major_00/';root=V/('all_lineage_Spectra_state_discovery_20261014_v2' if role=='Spectra' else 'all_lineage_cNMF_state_contrast_20261014_v2')
 for name,h in ledger['imported_file_sha256'][role].items():
  if name.startswith(prefix):assert sha(root/name)==h
 done.add((role,scope));record(role,scope,'REUSED_EXACT_COPY_PRIOR_NUMERIC_AUDIT',dict(producer_sha256=sha(producer),audit_sha256=sha(audit)))
def run(kind,scope):
 e=env.copy()
 if kind=='Spectra':
  diagnostic=['/tmp/diagnose_all_lineage_spectra_v2.py','--block',scope,'--plan',os.environ['PANCANCER_ALL_LINEAGE_SPECTRA_DIAGNOSTIC_PLAN']];audit='/tmp/audit_all_lineage_spectra_diagnostics_v2.py';e.update(PANCANCER_STATE_AUDIT_BLOCK=scope,PANCANCER_STATE_AUDIT_PLAN=os.environ['PANCANCER_ALL_LINEAGE_SPECTRA_DIAGNOSTIC_PLAN'])
 else:
  diagnostic=['/tmp/diagnose_all_lineage_cnmf_v2.py','--major',scope,'--plan',os.environ['PANCANCER_ALL_LINEAGE_CNMF_DIAGNOSTIC_PLAN']];audit='/tmp/audit_all_lineage_cnmf_diagnostics_v3.py';e['PANCANCER_STATE_AUDIT_MAJOR']=scope
 tag=kind+'_'+scope.replace(' ','_').replace('/','_');record(kind,scope,'DIAGNOSTIC_STARTED')
 with (O/(tag+'_diagnostic.log')).open('xb') as f:
  p=subprocess.Popen([os.sys.executable,*diagnostic],env=e,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT);record(kind,scope,'DIAGNOSTIC_PROCESS',dict(pid=p.pid));code=p.wait()
 if code:record(kind,scope,'DIAGNOSTIC_FAILED',dict(exit_code=code));raise RuntimeError(tag+' diagnostic failed; retained no retry')
 with (O/(tag+'_audit.log')).open('xb') as f:
  p=subprocess.Popen([os.sys.executable,audit],env=e,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT);record(kind,scope,'AUDIT_PROCESS',dict(pid=p.pid));code=p.wait()
 if code:record(kind,scope,'AUDIT_FAILED',dict(exit_code=code));raise RuntimeError(tag+' audit failed; retained no retry')
 done.add((kind,scope));record(kind,scope,'NUMERIC_DIAGNOSTIC_AUDIT_COMPLETED_NOT_STATE_PROMOTION')
while len(done)<24:
 progressed=False
 progress=json.loads((V/'all_lineage_Spectra_state_discovery_20261014_v2/progress.json').read_text());complete={(r['block'],r['seed']) for r in progress['summary']}
 for block in sp['blocks']:
  if ('Spectra',block) not in done and all((block,seed) in complete for seed in sp['seeds']):run('Spectra',block);progressed=True
 progress=json.loads((V/'all_lineage_cNMF_state_contrast_20261014_v2/progress.json').read_text());complete={r['major'] for r in progress['summary'] if r['all_consensus_completed']}
 for major in majors:
  if ('cNMF',major) not in done and major in complete:run('cNMF',major);progressed=True
 if len(done)==24:break
 alive=False;terminal=[]
 for job,producer in [('all_lineage_spectra_training_20261014_v2','train_all_lineage_spectra_v2.py'),('all_lineage_cnmf_training_20261014_v2','train_all_lineage_cnmf_v2.py')]:
  j=json.loads((V/(job+'.execution.json')).read_text());p=Path('/proc',str(j['child_pid']),'cmdline')
  if j['status']=='RUNNING':assert p.exists() and producer in p.read_bytes().decode(errors='replace');alive=True
  else:terminal.append(dict(job=job,status=j['status'],exit_code=j.get('exit_code')))
 if not alive:record('all','remaining','FITS_TERMINAL_WITH_INCOMPLETE_DIAGNOSTIC_SCOPE',terminal);raise RuntimeError('No live fit and incomplete scope; no restart')
 if not progressed:time.sleep(30)
receipt=dict(status='COMPLETE_ALL6_BLOCKS_AND18_MAJORS_FIXED_NUMERIC_DIAGNOSTICS_NOT_STATES',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),completed_scopes=len(done),clinical_LR_p_used=False,formal_states=False,original_labels_changed=False,remaining_goal=['fixed_programme_review_and_UCell_all_sources','all9CM_state_LR_receiver_response','state-awareDIALOGUE','state-aware_spatial']);(O/'inventory.json').write_text(json.dumps(receipt,indent=2)+chr(10));print(json.dumps(receipt),flush=True)
