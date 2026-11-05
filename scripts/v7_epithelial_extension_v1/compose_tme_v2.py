"""Explicit epithelial replacement overlay on concurrent TME v2; preserve Lodi panels."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

REPLACEMENTS = {
 'Epi_colonocyte_like':['Epi_colonocyte_like'],
 'Epi_goblet_like':['Epi_intestinal_goblet_like'],
 'Epi_intestinal_stem_like':['Epi_intestinal_stem_like'],
 'Epi_gastric_pit_mucous_like':['Epi_gastric_pit_like','Epi_gastric_gland_mucous_like'],
 'Epi_gastric_chief_parietal_like':['Epi_gastric_chief_like','Epi_gastric_parietal_like'],
 'Epi_enteroendocrine_like':['Epi_neuroendocrine_program'],
 'Epi_pancreatic_acinar_like':['Epi_pancreatic_acinar_like'],
 'Epi_pancreatic_ductal_like':['Epi_pancreatic_duct_like'],
 'Epi_islet_like':['Epi_pancreatic_islet_like'],
 'Epi_breast_luminal_like':['Epi_breast_LumHR_like','Epi_breast_LumSec_like'],
 'Epi_breast_basal_myoepithelial_like':['Epi_breast_basal_myo_like'],
 'Epi_prostate_basal_club_like':['Epi_prostate_basal_like','Epi_prostate_club_like'],
 'Epi_urothelial_like':['Epi_urothelial_differentiated_like'],
 'Epi_renal_distal_collecting_podocyte_like':['Epi_kidney_TAL_like','Epi_kidney_principal_like','Epi_kidney_podocyte_like'],
 'Epi_thyroid_follicular_like':['Epi_thyroid_follicular_like'],
 'Epi_mullerian_secretory_like':['Epi_fallopian_secretory_like'],
 'Epi_endometrial_like':['Epi_endometrial_secretory_like'],
 'Epi_squamous_suprabasal_like':['Epi_squamous_suprabasal_like'],
 'Epi_keratinocyte_like':['Epi_keratinocyte_differentiated_like'],
}

def compose(v2, epithelial):
    if len(v2['panels'])!=179: raise ValueError('Unexpected TME v2 scope; review overlay before applying')
    result=copy.deepcopy(v2)
    old={p['id']:p for p in v2['panels']}
    new={p['id']:p for p in epithelial['panels']}
    if not set(REPLACEMENTS)<=set(old): raise ValueError('TME v2 source changed; review mapping')
    for pid in REPLACEMENTS:
        if old[pid].get('source_layer')!='TME_de_novo_v1':
            raise ValueError('Refuse to replace author-code or unrelated source panel')
    added=list(dict.fromkeys(pid for ids in REPLACEMENTS.values() for pid in ids))
    if len(added)!=25 or not set(added)<=set(new): raise ValueError('Unexpected extension scope')
    result['panels']=[p for p in result['panels'] if p['id'] not in REPLACEMENTS]
    for p in result['panels']:
        if p['id']=='Epi_renal_tubular_like': p.update(copy.deepcopy(new[p['id']]))
    retained={p['id'] for p in result['panels']}
    if retained&set(added): raise ValueError('Unexpected retained panel ID collision')
    for pid in added:
        p=copy.deepcopy(new[pid]);p['source_layer']='V7_epithelial_literature_review_v1'
        p['source_locator']='resources/annotation/v7_epithelial_extension_v1/active_gene_evidence.json'
        result['panels'].append(p)
    refs={r['id']:r for r in result['references']}
    for ref in epithelial['references']:
        if ref['id'] not in refs:
            result['references'].append(copy.deepcopy(ref));refs[ref['id']]=result['references'][-1]
        elif refs[ref['id']]!=ref:
            if ref['id'] not in {'Kumar2023','Lake2023'}:
                raise ValueError('Unexpected reference conflict: '+ref['id'])
            # Same publications, different URLs/scope descriptions; preserve both.
            refs[ref['id']]['epithelial_review_reference']=copy.deepcopy(ref)
    result['policy'].update({k:v for k,v in epithelial['policy'].items()
                             if k in {'origin_assignment','benign_assignment','source27_claim'}})
    result['resource_version']='tme_v2_epithelial_review_v1'
    result['epithelial_overlay']={'replacements':REPLACEMENTS,'updated_retained_panel':'Epi_renal_tubular_like',
                                 'new_module_count':25,'retired_module_count':19,'final_cell_label_changes':0}
    return result

def write(repo,out):
    vp=repo/'resources/annotation/tme_annotation_panel_v2/registry.json'
    ep=repo/'resources/annotation/v7_epithelial_extension_v1/registry_merged.json'
    v2=json.loads(vp.read_text());epi=json.loads(ep.read_text());result=compose(v2,epi)
    out.mkdir(parents=True,exist_ok=False)
    rp=out/'tme_v2_registry_reviewed.json';rp.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(source_v2_sha256=sha(vp),source_epithelial_sha256=sha(ep),output_sha256=sha(rp),
                  source_v2_panels=len(v2['panels']),reviewed_panels=len(result['panels']),
                  replaced_panels=19,added_panels=25,updated_retained_panels=['Epi_renal_tubular_like'],
                  unchanged_panel_count=159,lodi_author_code_and_profile_panels_preserved=41,
                  original_files_modified=False,cell_label_changes=0,replacements=REPLACEMENTS)
    (out/'tme_v2_overlay_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();print(json.dumps(write(a.repo,a.out),indent=2))
