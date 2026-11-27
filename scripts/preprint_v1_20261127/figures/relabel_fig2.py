"""Replace internal version labels in the reused Figure 2 PDF with publication names (label text only; data marks untouched)."""
import pymupdf, re, sys
src, out = sys.argv[1], sys.argv[2]
doc = pymupdf.open(src); page = doc[0]
d = page.get_text('dict'); edits = []
for b in d['blocks']:
    for l in b.get('lines', []):
        text = ''.join(s['text'] for s in l['spans']); s0 = l['spans'][0]
        new = text
        new = re.sub(r'L-CM0?(\d)', lambda m: f'DM{int(m.group(1)):02d}', new)
        new = re.sub(r'R4 CM(\d+)', r'RM\1', new)
        new = new.replace('Legacy modules', 'Discovery modules').replace('R4 modules', 'Reference modules')
        if re.fullmatch(r'CM12', text.strip()): new = 'RM12'
        if new != text: edits.append((pymupdf.Rect(l['bbox']), new, s0['size'], s0['color'], s0['flags'] & 16, l['dir']))
for r, *_ in edits: page.add_redact_annot(r, fill=(1, 1, 1))
page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_NONE)
for r, new, size, color, bold, direc in edits:
    col = tuple(((color >> k) & 255) / 255 for k in (16, 8, 0))
    font = 'hebo' if bold else 'helv'
    w = pymupdf.get_text_length(new, fontname=font, fontsize=size)
    x = r.x0 + (r.width - w) / 2 if re.fullmatch(r'(RM|DM)\d+', new) else r.x0
    page.insert_text((x, r.y1 - size * 0.22), new, fontsize=size, fontname=font, color=col)
    print(f'{new!r} size={size:.1f}')
doc.save(out); print('saved', out, len(edits), 'labels')
