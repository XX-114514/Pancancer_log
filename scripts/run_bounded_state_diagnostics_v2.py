"""Explicit argument harness for immutable, scope-bound scientific diagnostics."""
import os,sys,runpy
mode=os.environ['PANCANCER_STATE_DIAGNOSTIC_MODE'];scope=os.environ['PANCANCER_STATE_DIAGNOSTIC_SCOPE'];plan=os.environ['PANCANCER_STATE_DIAGNOSTIC_PLAN']
if mode=='Spectra':script='/tmp/diagnose_all_lineage_spectra_v2.py';key='--block'
elif mode=='cNMF':script='/tmp/diagnose_all_lineage_cnmf_v2.py';key='--major'
else:raise ValueError(mode)
sys.argv=[script,key,scope,'--plan',plan];runpy.run_path(script,run_name='__main__')
