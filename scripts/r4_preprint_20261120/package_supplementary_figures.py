#!/usr/bin/env python3
"""Append externally supplied English legends to vector supplementary figures.

The private content file is never embedded in this generic packaging script.
Run with the Codex bundled Python runtime (pypdf and reportlab).
"""
import argparse
import io
import json
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph


def build(content, figures, target):
    writer = PdfWriter()
    for item in content['supplementary_figures']:
        source = PdfReader(str(figures / (item['id'] + '.pdf')))
        assert len(source.pages) == 1
        figure = source.pages[0]
        width, height = float(figure.mediabox.width), float(figure.mediabox.height)
        style = ParagraphStyle('Legend', fontName='Helvetica', fontSize=10,
                               leading=13, textColor='#000000')
        label = item['id'].split('_')[0].replace('FigureS', 'Figure S')
        text = str(item['en']).replace('\u2011', '-').replace('\u2013', '-')
        paragraph = Paragraph('<b>' + label + '.</b> ' + text, style)
        _, legend_height = paragraph.wrap(width - 48, 1000)
        band = legend_height + 36
        page = writer.add_blank_page(width=width, height=height + band)
        page.merge_transformed_page(figure, Transformation().translate(0, band))
        memory = io.BytesIO()
        overlay = canvas.Canvas(memory, pagesize=(width, height + band))
        paragraph.drawOn(overlay, 24, 18)
        overlay.save()
        memory.seek(0)
        page.merge_page(PdfReader(memory).pages[0])
    writer.add_metadata({'/Title': 'Supplementary figures for the R4 manuscript',
                         '/Author': ''})
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('wb') as handle:
        writer.write(handle)
    assert len(PdfReader(str(target)).pages) == len(content['supplementary_figures'])
    print(target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--content', type=Path, required=True)
    parser.add_argument('--figures', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(json.loads(args.content.read_text()), args.figures, args.output)
