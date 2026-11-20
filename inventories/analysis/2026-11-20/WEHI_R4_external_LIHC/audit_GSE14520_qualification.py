import csv,gzip,hashlib,json,sys,collections,datetime,pathlib,math
base,bulk,out=map(pathlib.Path,sys.argv[1:])
out.mkdir(parents=True,exist_ok=True)
m=base/'raw/cohorts/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz'
a=base/'metadata/platform_annotations/GPL3921.annot.gz'
s=bulk/'candidate_signatures.tsv'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(s)=='98907c4931b027282fde7a0a5a347b7ba36e5945f992ccad522555f31bddd5b2'
meta=collections.defaultdict(list); probes=set(); finite=True; n=0
with gzip.open(m,'rt') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row: continue
  if row[0].startswith('!Sample_'): meta[row[0]].append(row[1:])
  elif row[0]=='ID_REF': n=len(row)-1
  elif not row[0].startswith('!') and n:
   assert len(row)==n+1 and row[0] not in probes
   probes.add(row[0]); finite &= all(math.isfinite(float(v)) for v in row[1:])
assert finite and n==len(meta['!Sample_geo_accession'][0])
fields=collections.Counter(); tissue=collections.Counter(); individuals=[]
for row in meta['!Sample_characteristics_ch1']:
 assert len(row)==n
 for v in row:
  k,sep,val=v.partition(':'); k=k.strip().lower(); fields[k]+=1
  if k=='tissue': tissue[val.strip().lower()]+=1
  if k=='individual': individuals.append(val.strip())
unique=collections.defaultdict(set); ambiguous=collections.defaultdict(set)
with gzip.open(a,'rt') as f:
 for line in f:
  if line.startswith('ID\t'): header=line.rstrip('\n').split('\t'); break
 for row in csv.DictReader(f,fieldnames=header,delimiter='\t'):
  if row['ID'].startswith('!'): break
  if row['ID'] not in probes: continue
  genes={g.strip() for g in row['Gene symbol'].split('///') if g.strip() and g.strip()!='---'}
  for g in genes: (unique if len(genes)==1 else ambiguous)[g].add(row['ID'])
rows=[]
with s.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['cm']=='CM13': rows.append(dict(branch=r['branch'],cm=r['cm'],gene=r['gene'],unambiguous_probes=len(unique[r['gene']]),ambiguous_probes=len(ambiguous[r['gene']])))
with (out/'CM13_probe_coverage.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
summary={b:dict(planned=sum(r['branch']==b for r in rows),unambiguous=sum(r['branch']==b and r['unambiguous_probes']>0 for r in rows),missing=[r['gene'] for r in rows if r['branch']==b and not r['unambiguous_probes']]) for b in sorted({r['branch'] for r in rows})}
report=dict(status='EXPRESSION_AUDITED_CLINICAL_NOT_QUALIFIED',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in (m,a,s)},samples=n,probe_rows=len(probes),all_values_finite=finite,characteristic_fields=dict(fields),tissue_counts=dict(tissue),unique_individual_labels=len(set(individuals)),individual_labels_are_not_verified_clinical_bridge=True,coverage=summary,clinical_endpoint_available_in_matrix=False,clinical_patient_join_verified=False,survival_model_run=False,previously_explored_dataset=True,annotation_policy='Exact official gene symbol; multi-symbol probes excluded; no alias rescue or outcome-based probe choice',output_sha256=sha(out/'CM13_probe_coverage.tsv'))
(out/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
