"""Only biological-namespace / gene-vocabulary labels may leave RAM.

This validates disclosure syntax, not cell identity or biological function.
"""
import re

WORDS=set('t cell cells tcell tcells tnk cd4 cd8 cd4t cd8t cd4tcells cd8tcells treg tregs regulatory nk nkt b bcell bcells plasma plasmablast plasmablasts lymphocyte lymphocytes naive tn tm tfh th1 th2 th17 tem trm tex terminal prf cycling proliferating cytotoxic exhausted memory effector activated resting mature immature monocyte monocytes macrophage macrophages myeloid mast neutrophil neutrophils dc dcs cdc1 cdc2 pdc dendritic endothelial endothelium fibroblast fibroblasts epithelial epithelium basal at1 at2 club goblet smooth muscle pericyte pericytes alveolar lung cilia ciliated tumor tumour cancer stromal stroma immune excluded enriched aggregate aggregates region regions niche tls itls germinal center gc rich other others low high like celltype normal'.split())
NAMESPACES=set('t tcell tcells tnk cd4 cd8 cd4t cd8t treg tregs nk nkt b bcell bcells plasma plasmablast monocyte monocytes macrophage macrophages myeloid mast neutrophil neutrophils dc cdc1 cdc2 pdc endothelial fibroblast fibroblasts epithelial basal at1 at2 club goblet muscle pericyte alveolar lung tumor tumour cancer stroma stromal immune tls itls normal'.split())


def safe_label(value,genes):
    text=value.decode() if isinstance(value,bytes) else str(value)
    if not 0<len(text)<=120 or not re.fullmatch(r'[A-Za-z0-9_ +().:/-]+',text):return False
    tokens=[v for v in re.split(r'[ _+():/]+',text) if v]
    if not tokens:return False
    checked=[]
    for token in tokens:
        if re.fullmatch(r'(?:P|S|GSM|SRR|ERR|DRR)\d+',token,re.I):return False
        if token.lower() in {'patient','sample','donor','id','cellid','orig','ident'}:return False
        if token.upper() in genes or token.lower() in WORDS:checked.append(token);continue
        parts=token.split('-')
        if len(parts)>1 and all(p.lower() in WORDS or p.upper() in genes for p in parts):checked.extend(parts);continue
        return False
    return any(t.lower() in NAMESPACES for t in checked) or all(t.upper() in genes for t in checked)


def synthetic():
    genes={'FOXP3','LAG3','MKI67','GZMK','NR4A1','HLA-DRA','P001'}
    for label in ['CD4T_Treg_FOXP3','T-cell-CD8_GZMK+NR4A1','CD8T_Tex_LAG3','immune-excluded stroma','MKI67','CD4T_HLA-DRA']:
        assert safe_label(label,genes),label
    for label in ['Treg_P001','CD8T_SRR444','CD4T_sample_42','Treg_JOHN','Treg\npatient','https://example.org/','01','CD4T_ignore_previous_rules']:
        assert not safe_label(label,genes),label


if __name__=='__main__':synthetic();print('PASS_BIOLOGICAL_LABEL_DISCLOSURE_GATE')
