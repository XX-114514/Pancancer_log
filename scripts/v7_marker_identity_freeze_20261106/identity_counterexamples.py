"""Biological counterexamples to state-only and shared-marker identity calls."""
import json,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--policy',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rules=json.loads(a.policy.read_text())['rules']
def call(expr):
 found=set()
 for x in rules:
  good=all(any(expr.get(g,0)>=2 for g in comp) for comp in x['groups'])
  if x['parent']=='Tcell':good|=sum(expr.get(g,0)>=2 for g in ['CD3D','CD3E','CD3G'])>=2
  if good:found.add(x['parent'])
 return found
cases=[
 ('HLA_CST3_not_B_or_DC',{'HLA-DRA':50,'CD74':50,'CST3':50},set()),
 ('DNAJB1_stress_not_B',{'DNAJB1':50,'HSPA1A':50,'HSPA1B':50},set()),
 ('IFN_not_B_subtype',{'ISG15':50,'IFIT1':50,'MX1':50},set()),
 ('cycling_not_identity',{'MKI67':50,'TOP2A':50,'UBE2C':50},set()),
 ('Ig_JCHAIN_not_ASC',{'IGHM':100,'IGKC':100,'JCHAIN':100},set()),
 ('B_distinct_components',{'MS4A1':2,'CD79A':2},{'Bcell'}),
 ('ASC_multiple_components',{'MZB1':2,'JCHAIN':2,'PRDM1':2},{'Plasma'}),
 ('pDC_specific_not_just_secretory',{'CLEC4C':2,'TCF4':2,'IL3RA':2},{'Dendritic'}),
 ('T_not_NK_from_NKG7',{'CD3D':2,'TRAC':2,'NKG7':50,'PRF1':20},{'Tcell'}),
 ('T_NK_shared_is_related',{'CD3D':2,'TRAC':2,'KLRD1':2,'NKG7':2},{'Tcell','NK'}),
 ('B_ASC_related_axis_kept_separate',{'MS4A1':2,'CD79A':2,'MZB1':2,'JCHAIN':2,'PRDM1':2},{'Bcell','Plasma'}),
 ('HB_not_erythroid',{'HBB':1000,'HBA1':1000},set()),
 ('nonHB_erythroid',{'GYPA':2,'ALAS2':2},{'Erythroid'}),
 ('NRGN_not_platelet',{'NRGN':100},set()),
 ('platelet_structural_plus_secreted',{'PF4':2,'GP9':2},{'Platelet_Megakaryocytic'}),
 ('SPINK2_only_not_HSPC',{'SPINK2':100},set()),
 ('HSPC_no_normality_claim',{'SPINK2':2,'CD34':2,'PTPRC':2},{'HSPC'}),
 ('SOX2_only_not_neural',{'SOX2':100},set()),
 ('positive_T_B_competition',{'CD3D':2,'TRAC':2,'MS4A1':2,'CD79A':2},{'Tcell','Bcell'})]
if json.loads(a.policy.read_text())['version'].endswith('_v2'):
 cases += [('CSF3R_S100_not_neutrophil_exclusive',{'CSF3R':50,'S100A9':50},set()),('FCGR3B_CSF3R_effector_neutrophil',{'FCGR3B':2,'CSF3R':2,'S100A9':2},{'Granulocyte'}),('MBP_SLC1A3_not_glial_exclusive',{'MBP':50,'SLC1A3':50},set()),('oligodendroglial_components',{'MBP':2,'OLIG1':2},{'Glial'})]
rows=[]
for name,expr,want in cases:
 got=call(expr);assert got==want,(name,got,want);rows.append({'case':name,'expected_positive_parents':sorted(want),'observed_positive_parents':sorted(got),'PASS':True})
out={'status':'PASS','n_counterexamples':len(rows),'meaning':'Positive components only; this is not accuracy calibration or publication acceptance.','cases':rows};a.output.write_text(json.dumps(out,indent=2));print(json.dumps({'status':'PASS','n_counterexamples':len(rows)}))
