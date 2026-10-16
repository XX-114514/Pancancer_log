"""Conservative source-author text roles, never a replacement annotation."""
import re

VIEWS=['author_CD8','author_Treg','author_Macrophage','author_B_Plasma']
CHANNELS=['frozen_candidates','author_unambiguous','author_concordant']
PARENTS={
    'author_CD8':{'T','T cell','T/NK lymphoid','Lymphoid cell'},
    'author_Treg':{'T','T cell','T/NK lymphoid','Lymphoid cell'},
    'author_Macrophage':{'Macrophage','Myeloid'},
    'author_B_Plasma':{'B cell','B','Plasma','Lymphoid cell'},
}


def author_role(text):
    s=re.sub(r'[_\-]+',' ',str(text).strip().casefold())
    s=re.sub(r'\s+',' ',s)
    if s in {'','nan','none','unknown','not reviewed','unavailable','unassigned'}:return -1
    if re.search(r'\bb\s+like\s+t\s*cells?\b',s):return -2
    hit=[]
    if re.search(r'\bcd8(?:\b|t)',s):hit.append(0)
    if re.search(r'\btregs?\b|\bregulatory\s+t(?:\b|cells?)|\bt\s+regulatory\b',s):hit.append(1)
    if re.search(r'\bmacrophages?\b|\bmacro\b|\btams?\b|mφ',s):hit.append(2)
    if re.search(r'\bb\s*cells?\b|\bb\s+lymph|\bplasma(?:blasts?|\s+cells?)?\b|\bantibody\s+secreting\s+b\b',s) or s in {'b','b lymphocyte','b lymphocytes'}:hit.append(3)
    if len(set(hit))>1:return -2
    if not hit:return -1
    role=hit[0]
    if role==0 and re.search(r'\bcd4(?:\b|t)',s):return -2
    if role==2 and re.search(r'\bmonocytes?\b|\bmono\b|\bdendritic\b|\bdc\b|\bcdc[12]?\b',s):return -2
    if role==3 and re.search(r'\bt\s*cells?\b|\bcd[48](?:\b|t)',s):return -2
    return role


def phase(text):
    s=str(text).strip().casefold().replace('-','_').replace(' ','_')
    if s in {'pre','pre_treatment','pretreatment','baseline','pre_vaccination'}:return 'pre'
    if s in {'post','post_treatment','posttreatment','post_vaccination'}:return 'post'
    if s in {'on_trt','on_treatment','ontreatment'}:return 'on_treatment'
    if s in {'on_prog','progression','post2','recurrence'}:return 'progression_or_recurrence'
    if s in {'','unknown','not_reviewed','nan','none'}:return 'unknown'
    return 'other_recorded_phase'


def material(row):
    tissue=str(row.tissue).strip().casefold();sample=str(row.sample_class).strip().casefold();capture=str(row.capture_scope).strip().casefold()
    if 'conflict' in tissue or 'conflict' in sample:return 'CONFLICTING_METADATA'
    if str(row.ex_vivo).casefold()=='true':return 'EX_VIVO'
    if 'infusion' in tissue or 'infusion' in capture or 'car_t' in capture:return 'INFUSION_PRODUCT'
    if 'blood' in tissue or tissue in {'pbmc','pbmcs','peripheral blood'}:return 'BLOOD'
    if 'marrow' in tissue:return 'BONE_MARROW'
    if tissue in {'ln','lymph_node','lymph node'} or 'lymph node' in tissue or 'lymph_node' in tissue:return 'LYMPH_NODE'
    if 'normal' in tissue or 'adjacent' in tissue or 'normal' in sample or 'adjacent' in sample:return 'NORMAL_OR_ADJACENT'
    if 'metasta' in tissue or 'metasta' in sample:return 'RECORDED_SOLID_METASTASIS'
    if 'primary tumor' in tissue or 'primary tumour' in tissue or ('primary' in sample and ('tumor' in sample or 'tumour' in sample)):return 'RECORDED_PRIMARY_TUMOR'
    if any(x in tissue+' '+sample for x in ['tumor','tumour','cancer']):return 'RECORDED_SOLID_TUMOR'
    if 'whole_tumor' in capture or 'serial_tumor_biopsies' in capture:return 'TUMOR_CAPTURE_SITE_UNSPECIFIED'
    if tissue not in {'','unknown','not_reviewed','none','nan'}:return 'ANATOMICAL_SITE_ONLY'
    return 'UNKNOWN_MATERIAL'


def rule_fixture():
    expected={'CD8+ T cells':0,'CD8T_ISG15':0,'CD4/CD8 T cells':-2,'Treg':1,'regulatory T cell':1,
              'macrophage':2,'Macro_APOE':2,'monocyte/macrophage':-2,'macrophage/DC':-2,
              'B cells':3,'Bcell':3,'plasma cells':3,'Antibody_secreting_B_cell':3,'B-like T cells':-2,
              'T/NK':-1,'Naive T cells':-1,'Fibroblast':-1,'unknown':-1}
    for name,value in expected.items():assert author_role(name)==value,(name,value)
