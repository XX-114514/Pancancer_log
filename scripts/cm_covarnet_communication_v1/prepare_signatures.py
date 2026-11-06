from pathlib import Path
import csv,json,zipfile,xml.etree.ElementTree as E,collections,gzip,hashlib
p=Path(__file__).resolve().parents[1];h2m=collections.defaultdict(set);m2h=collections.defaultdict(set)
for row in csv.reader((p/'sources/HMD_HumanPhenotype.rpt').open(),delimiter='\t'):
 if len(row)>=3 and row[0] and row[2]:h2m[row[0]].add(row[2]);m2h[row[2]].add(row[0])
mapping={m:next(iter(hs)) for m,hs in m2h.items() if len(hs)==1 and len(h2m[next(iter(hs))])==1}
z=zipfile.ZipFile(p/'sources/immune_dictionary_S3.xlsx');ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};ss=[''.join(e.itertext()) for e in E.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)];sheets=E.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns)
rows=[];unmapped=collections.Counter();sign=collections.Counter()
for i,s in enumerate(sheets):
 data=E.fromstring(z.read(f'xl/worksheets/sheet{i+1}.xml')).find('m:sheetData',ns)
 for row in list(data)[1:]:
  vals=[ss[int(c.find('m:v',ns).text)] if c.get('t')=='s' else (c.find('m:v',ns).text if c.find('m:v',ns) is not None else '') for c in row]
  ct,_,cyt,cytstr,mg,fc,fdr,*_=vals;direction='up' if float(fc)>0 else 'down';sign[direction]+=1
  if mg not in mapping:unmapped[ct]+=1;continue
  rows.append([ct.replace(' ','_'),s.attrib['name'],cyt,cytstr,mg,mapping[mg],float(fc),float(fdr),direction])
with gzip.open(p/'sources/signatures_one_to_one.tsv.gz','wt') as handle:
 w=csv.writer(handle,delimiter='\t');w.writerow(['cell_type','sheet','cytokine','cytokine_symbol_mouse','mouse_gene','human_gene','log2FC','FDR','direction']);w.writerows(rows)
(p/'config/signature_audit.json').write_text(json.dumps({'reciprocal_one_to_one_mappings':len(mapping),'mapped_signature_rows':len(rows),'original_directions':dict(sign),'unmapped_by_type':dict(unmapped),'orthology_source':'current_MGI_20261106_not_author_archived_snapshot','author_direction_handling':'code selects Cytokine,Gene only; project preserves signs'},indent=2))
print((p/'config/signature_audit.json').read_text())
