"""Minimal Markdown -> DOCX for the preprint draft (headings, paragraphs, bullets, numbered lists, tables, bold/italic, figures)."""
import re, sys
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
src, out, figdir = sys.argv[1], sys.argv[2], sys.argv[3]
FIGS = {'Figure 1.': 'Fig1_release_lineage.png', 'Figure 2.': 'Fig2_source_dependence_reused_v12.png', 'Figure 3.': 'Fig3_module_discovery_audit.png',
        'Figure 4.': 'Fig4_R4_module_catalogue.png', 'Figure 5.': 'Fig5_vascular_positive_control_reused_v12.png', 'Figure 6.': 'Fig6_Tex_centred_module.png',
        'Supplementary Figure S1.': 'FigS_bulk_readout_v2.png'}
doc = Document(); st = doc.styles['Normal']; st.font.name = 'Arial'; st.font.size = Pt(10)
def runs(p, text):
    for part in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)', text):
        if not part: continue
        if part.startswith('**'): p.add_run(part[2:-2]).bold = True
        elif part.startswith('*') and part.endswith('*'): p.add_run(part[1:-1]).italic = True
        elif part.startswith('^'): r = p.add_run(part[1:-1]); r.font.superscript = True
        else: p.add_run(part)
lines = open(src, encoding='utf-8').read().split('\n'); i = 0
if lines[0] == '---':
    j = lines.index('---', 1); meta = {l.split(':', 1)[0]: l.split(':', 1)[1].strip().strip('"') for l in lines[1:j] if ':' in l}
    h = doc.add_heading(meta.get('title', ''), 0); p = doc.add_paragraph(); p.add_run(meta.get('subtitle', '')).italic = True; i = j + 1
pending_fig = None
while i < len(lines):
    l = lines[i]
    if l.startswith('|'):
        rows = []
        while i < len(lines) and lines[i].startswith('|'):
            if not re.match(r'^\|[\s\-|:]+\|$', lines[i]): rows.append([c.strip() for c in lines[i].strip('|').split('|')])
            i += 1
        t = doc.add_table(rows=len(rows), cols=len(rows[0])); t.style = 'Table Grid'
        for r, row in enumerate(rows):
            for c, val in enumerate(row):
                cell = t.cell(r, c); cell.text = ''; pp = cell.paragraphs[0]; runs(pp, val)
                for rr in pp.runs: rr.font.size = Pt(8); rr.bold = rr.bold or r == 0
        continue
    if l.startswith('#'):
        lvl = len(l) - len(l.lstrip('#')); doc.add_heading(l.lstrip('#').strip(), min(lvl, 3))
    elif l.strip() == '---' or not l.strip(): pass
    elif re.match(r'^\s*- ', l):
        ind = (len(l) - len(l.lstrip())) // 2; p = doc.add_paragraph(style='List Bullet' if ind == 0 else 'List Bullet 2'); runs(p, l.strip()[2:])
    elif re.match(r'^\s*\d+\. ', l):
        p = doc.add_paragraph(style='List Number'); runs(p, re.sub(r'^\s*\d+\. ', '', l))
    else:
        fig = next((v for k, v in FIGS.items() if l.startswith('**' + k)), None)
        if fig:
            doc.add_picture(f'{figdir}/{fig}', width=Inches(6.5)); doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        p = doc.add_paragraph(); runs(p, l.strip())
    i += 1
doc.save(out); print('saved', out)
