#!/usr/bin/env python3
"""Rebuild every manuscript figure using only the delivered plotting package."""
from pathlib import Path
import argparse,subprocess,sys,json,hashlib,shutil

def run(root,out):
 out.mkdir(parents=True,exist_ok=False);code=Path(__file__).resolve().parent;inputs=root/'inputs';work=out/'intermediate';work.mkdir()
 jobs=[('r4',[str(code/'render_R4_portable.py'),'--root',str(inputs/'r4'),'--run',str(inputs/'r4/R4_preprint_20261120_v1'),'--output',str(work/'r4')]),('composition',[str(code/'render_composition_and_LN.py'),'--root',str(inputs),'--prior',str(inputs/'R4_figure_tables'),'--configuration',str(inputs/'figure_configuration.json'),'--output',str(work/'composition')]),('story',[str(code/'render_story_panels.py'),'--inputs',str(inputs),'--output',str(work/'story')])]
 for name,args in jobs:
  with (out/(name+'.log')).open('w') as f:subprocess.run([sys.executable,'-B']+args,check=True,stdout=f,stderr=subprocess.STDOUT)
  print('COMPLETE',name,flush=True)
 mapping={
 'Figure1_resolved_atlas':('r4','Figure2_unified_atlas'),
 'Figure2_composition_landscape':('composition','Figure3_composition_landscapes'),
 'Figure3_immune_state_combination':('story','Figure3_immune_state_combination'),
 'Figure4_immune_clinical_readouts':('story','Figure4_immune_clinical_readouts'),
 'Figure5_lymph_node_organization':('composition','Figure5_T_state_spatial_context'),
 'Figure6_vascular_architecture':('r4','Figure4_spatial_organization'),
 'Figure7_cancer_context':('story','Figure7_cancer_context'),
 'FigureS01_source_coverage':('r4','FigureS01_complete_source_coverage'),
 'FigureS02_annotation_dotplot':('r4','FigureS02_annotation_dotplot'),
 'FigureS03_module_networks':('r4','FigureS03_all_module_networks'),
 'FigureS04_source_cancer_coordinates':('r4','FigureS04_source_and_cancer_embedding'),
 'FigureS05_spatial_denominators':('r4','FigureS05_spatial_denominator_sensitivity'),
 'FigureS06_identity_key':('r4','FigureS10_full_identity_UMAP'),
 'FigureS07_immune_atlas':('story','FigureS07_immune_atlas'),
 'FigureS08_T_state_dimensions':('story','FigureS08_T_state_dimensions'),
 'FigureS09_prediction_information':('story','FigureS09_prediction_information'),
 'FigureS10_source_composition':('composition','FigureS14_legacy_source_contexts'),
 'FigureS11_imaging_neighborhoods':('r4','Figure5_imaging_neighborhoods'),
 'FigureS12_molecular_context':('composition','Figure8_communication_across_designs')}
 target=out/'figures';target.mkdir();receipt=[]
 for name,(branch,source) in mapping.items():
  for ext in ['png','pdf','svg']:
   p=work/branch/(source+'.'+ext);q=target/(name+'.'+ext);shutil.copy2(p,q)
   receipt.append({'file':q.name,'producer':branch,'intermediate':p.name,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
 (out/'figure_manifest.json').write_text(json.dumps({'status':'ALL_MANUSCRIPT_FIGURES_REBUILT','figures':len(mapping),'files':receipt,'new_inferential_tests':False},indent=2)+'\n');print('ALL_FIGURES_REBUILT',len(mapping),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--package-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.package_root,a.output)
