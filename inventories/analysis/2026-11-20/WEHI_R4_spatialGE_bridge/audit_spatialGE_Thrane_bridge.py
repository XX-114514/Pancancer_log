import pathlib,sys,zipfile,hashlib,csv,io,json,datetime,collections
root=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
p=root/'spatialGE_Thrane_v1/melanoma_thrane.zip';part=p.with_suffix('.zip.part')
if not p.exists(): p=part
sha=lambda b:hashlib.sha256(b).hexdigest()
assert p.stat().st_size==3400344 and sha(p.read_bytes())=='c7957972360dd4753d6d49764d98d8a7f63c070ec2de4e4d968f2fc51881f32c'
z=zipfile.ZipFile(p);assert z.testzip() is None
original=root/'Thrane2018_histology_v1/Thrane_et_al_2018_CAN_RES.zip'
assert sha(original.read_bytes())=='8adc48e90b50c196fb9b8efffcd9cf8d72bc005bcdfcc50ab806a00627bb2c1e'
a=zipfile.ZipFile(original)
read=lambda archive,name:list(csv.reader(io.StringIO(archive.read(name).decode()),delimiter='\t'))
results=[]
for name in sorted(z.namelist()):
 if not(name.startswith('melanoma_thrane/') and name.endswith('_counts.tsv')):continue
 sample=pathlib.Path(name).name.removesuffix('_counts.tsv')
 oldname=next(n for n in a.namelist() if n.startswith('Thrane_et_al_2018_CAN_RES/Count_data/') and n.endswith(sample+'_counts.tsv'))
 old=read(a,oldname);new=read(z,name)
 assert old[0]==new[0] and len(old)==len(new)
 assert len(set(old[0][1:]))==len(old[0])-1
 mismatches=sum(x[1:]!=y[1:] for x,y in zip(old[1:],new[1:]))
 symbols=sum(x[0].split()[0]!=y[0] for x,y in zip(old[1:],new[1:]))
 assert mismatches==symbols==0
 mapping_name=name.replace('_counts.tsv','_mapping.tsv');mapping=list(csv.DictReader(io.StringIO(z.read(mapping_name).decode()),delimiter='\t'))
 assert {r['spot'] for r in mapping}==set(new[0][1:]) and len(mapping)==len(new[0])-1
 assert all(tuple(map(int,r['spot'].split('x')))==(int(r['xpos']),int(r['ypos'])) for r in mapping)
 results.append(dict(sample=sample,spots=len(mapping),genes=len(new)-1,numeric_values=(len(new)-1)*len(mapping),all_count_rows_exact=True,all_symbols_match_after_ensembl_suffix_removal=True,coordinates_exact_barcode_integer_grid=True,duplicate_symbols=len(new)-1-len({r[0] for r in new[1:]}),original_count_sha256=sha(a.read(oldname)),derivative_count_sha256=sha(z.read(name)),mapping_sha256=sha(z.read(mapping_name))))
clinical=list(csv.DictReader(io.StringIO(z.read('melanoma_thrane/thrane_clinical.csv').decode())))
assert {r['sampleID'] for r in clinical}=={r['sample'] for r in results}
assert len(clinical)==len(results)
report=dict(status='DERIVATIVE_BRIDGE_AUDITED_NOT_PRIMARY_PATIENT_VERIFIED',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),archive_sha256=sha(p.read_bytes()),archive_bytes=p.stat().st_size,zip_crc='PASS',samples=results,clinical_rows=len(clinical),distinct_derivative_patient_labels=len({r['patient'] for r in clinical}),clinical_fields=list(clinical[0]),clinical_member_sha256=sha(z.read('melanoma_thrane/thrane_clinical.csv')),primary_patient_bridge_verified=False,image_registration_provided=False,includes_images=False,clinical_outcomes_used=False,other_four_sections_patient_bridge=False,source='https://github.com/FridleyLab/spatialGE_Data',scope='All four packaged slices, all counts and coordinates; no extrapolation to absent slices')
(out/'complete.json').write_text(json.dumps(report,indent=2)+'\n')
if p==part:p.rename(p.with_suffix(''))
print(json.dumps(report,indent=2))
