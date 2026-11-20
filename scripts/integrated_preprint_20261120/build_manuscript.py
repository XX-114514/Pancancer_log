#!/usr/bin/env python3
"""Build editable bilingual manuscripts from an external private content JSON.

Use the Codex bundled Python runtime. Scientific prose and private figure data
are supplied separately and are intentionally not embedded in this code.
"""
import argparse,json,math
from pathlib import Path
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from PIL import Image


def text(s):
    return str(s).replace('\u2011','-').replace('\u2013','-').replace('\u2014','-').replace('\u2212','-')


def font(style,lang,size,bold=False):
    style.font.name='Times New Roman';style.font.size=Pt(size);style.font.bold=bold;style.font.color.rgb=RGBColor(0,0,0)
    rpr=style.element.get_or_add_rPr();rf=rpr.rFonts
    if rf is None:rf=OxmlElement('w:rFonts');rpr.insert(0,rf)
    rf.set(qn('w:ascii'),'Times New Roman');rf.set(qn('w:hAnsi'),'Times New Roman');rf.set(qn('w:eastAsia'),'Arial Unicode MS' if lang=='zh' else 'Times New Roman')
    for key in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme']:
        rf.attrib.pop(qn('w:'+key),None)
    ppr=style.element.pPr
    if ppr is not None:
        for item in list(ppr):
            if item.tag==qn('w:pBdr'):ppr.remove(item)


def build(content,figures,output):
    output.mkdir(parents=True,exist_ok=True)
    for lang in ['en','zh']:
        d=Document()
        for style in d.styles:
            ppr=style.element.find(qn('w:pPr'))
            if ppr is not None:
                for el in list(ppr):
                    if el.tag==qn('w:pBdr'):ppr.remove(el)
        sec=d.sections[0]
        sec.page_width=Inches(8.5);sec.page_height=Inches(11)
        sec.top_margin=sec.bottom_margin=Inches(.85)
        sec.left_margin=sec.right_margin=Inches(.85)
        sec.header_distance=sec.footer_distance=Inches(.35)
        font(d.styles['Normal'],lang,10.5 if lang=='zh' else 11);s=d.styles['Normal'].paragraph_format;s.line_spacing=1.1 if lang=='zh' else 1.15;s.space_after=Pt(5 if lang=='zh' else 6);s.widow_control=True
        font(d.styles['Title'],lang,18,True);d.styles['Title'].paragraph_format.space_after=Pt(14)
        for name,size in [('Heading 1',13),('Heading 2',11.5)]:
            font(d.styles[name],lang,size,True);d.styles[name].paragraph_format.keep_with_next=True;d.styles[name].paragraph_format.space_before=Pt(12);d.styles[name].paragraph_format.space_after=Pt(6)
        font(d.styles['Caption'],lang,9.5);d.styles['Caption'].paragraph_format.space_after=Pt(5);d.styles['Caption'].paragraph_format.line_spacing=1.0
        d.core_properties.author='';d.core_properties.last_modified_by='';d.core_properties.title=content['title_'+lang]
        d.add_paragraph(content['title_'+lang],'Title')
        for part in content['sections']:
            d.add_heading(text(part['title_'+lang]),level=part['level'])
            for para in part['paragraphs']:d.add_paragraph(text(para[lang]))
        if lang=='en':d.add_page_break()
        d.add_heading('References' if lang=='en' else '参考文献',1)
        for i,reference in enumerate(content['references'],1):
            p=d.add_paragraph(f'[{i}] '+text(reference));p.paragraph_format.left_indent=Inches(.22);p.paragraph_format.first_line_indent=Inches(-.22)
            p.paragraph_format.space_after=Pt(4)
            for run in p.runs:run.font.size=Pt(10)
        # Figures have complete legends and their own page; no raw data tables.
        for item in content['main_figures']:
            d.add_page_break();label=item['id'].split('_')[0].replace('Figure','Figure ' if lang=='en' else '图')
            p=d.add_paragraph(label,'Heading 1');p.paragraph_format.space_before=Pt(0);p.paragraph_format.space_after=Pt(6)
            f=figures/(item['id']+'.png')
            with Image.open(f) as im:w,h=im.size
            legend=text(item[lang]);lines=math.ceil(len(legend)/(105 if lang=='en' else 46))
            caption_height=max(.65,lines*12/72+.12)
            maxheight=min(7.1,8.45-caption_height)
            width=min(6.7,maxheight*w/h)
            p=d.add_paragraph();p.paragraph_format.space_after=Pt(6);p.paragraph_format.keep_with_next=True;p.alignment=1
            run=p.add_run();run.add_picture(str(f),width=Inches(width))
            pic=run._r.xpath('.//wp:docPr')
            if pic:pic[0].set('descr',legend)
            p=d.add_paragraph(legend,'Caption');p.paragraph_format.widow_control=True;p.paragraph_format.keep_together=True
        d.add_page_break();d.add_heading('Supplementary figure legends' if lang=='en' else '补充图图注',1)
        d.add_paragraph('Supplementary figures are supplied as PDFs and high-resolution images. Archived RNA-program panels retain their raster provenance.' if lang=='en' else '补充图以PDF和高分辨率图片保存在图包中；存档RNA程序图保留其栅格来源。科研图使用英文。')
        for item in content['supplementary_figures']:
            p=d.add_paragraph(item['id'].split('_')[0].replace('FigureS','Figure S' if lang=='en' else '图S'),'Heading 2');p.paragraph_format.keep_with_next=True;p.paragraph_format.space_before=Pt(6);p.paragraph_format.space_after=Pt(4)
            p=d.add_paragraph(text(item[lang]));p.paragraph_format.line_spacing=1.0;p.paragraph_format.space_after=Pt(4)
            for run in p.runs:run.font.size=Pt(10)
        footer=sec.footer.paragraphs[0];footer.alignment=1
        field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)
        target=output/(content.get('output_prefix','PanCancer_integrated_manuscript')+'_'+lang+'.docx');d.save(target)
        print(target,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--content',type=Path,required=True);p.add_argument('--figures',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();build(json.loads(a.content.read_text()),a.figures,a.output)
