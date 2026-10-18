"""Fetch only fixed, small author-public compartment inputs; no effects/models."""
import argparse, csv, hashlib, io, json, urllib.request
from pathlib import Path
from datetime import datetime, timezone

REVISION='7cb78392967696fc5a6ec44d05e687a11052f94f'
FILES={
 'CK_av_3B.csv':(5008370,'94dd626cb5ed9aa0c54ecbcaa1b6cbf5c4e6f211'),
 'Stroma_av_3B.csv':(3132965,'74e4f604b48eaa99e68d450eb6d0776d17737a58'),
 'UQ_CK_large_PFS_OS_only_1st_2nd_line.csv':(666536,'df16a07832edd52a766b5cb4424f86f42965ed40'),
 'UQ_Stroma_large_PFS_OS_only_1st_2nd_line.csv':(600703,'6002')}
# Exact Git blob from the source tree; no GitHub credentials or author model code.
FILES['UQ_Stroma_large_PFS_OS_only_1st_2nd_line.csv']=(600703,'a91b62a3dc3827ba8a5499f4563b6ab9b92b362a')

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 plan=json.loads(a.plan.read_text());genes=set(sum(plan['panels'].values(),[]));rows=[]
 for name,(size,blob) in FILES.items():
  url=f'https://raw.githubusercontent.com/tznaung/NSCLC_SpatialOmics/{REVISION}/data/{name}'
  req=urllib.request.Request(url,headers={'User-Agent':'CM02-public-source-qualification/1'})
  with urllib.request.urlopen(req,timeout=45) as response:data=response.read(16*1024**2)
  assert len(data)==size and len(data)<16*1024**2
  assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==blob
  (a.output/name).write_bytes(data)
  records=list(csv.reader(io.StringIO(data.decode('utf-8-sig'))));width=len(records[0])
  assert all(len(x)==width for x in records)
  head=records[0];first=[x[0] for x in records[1:]]
  col_hits=genes&set(head);row_hits=genes&set(first)
  orientation='genes_in_columns' if len(col_hits)>len(row_hits) else 'genes_in_rows'
  measured=col_hits if orientation=='genes_in_columns' else row_hits
  ambiguous=[g for g in genes if (head if orientation=='genes_in_columns' else first).count(g)>1]
  measured-=set(ambiguous)
  coverage={k:dict(required=len(gg),unique_measured=len(set(gg)&measured),missing=sorted(set(gg)-measured)) for k,gg in plan['panels'].items()}
  rows.append(dict(file=name,url=url,bytes=size,git_blob_sha1=blob,sha256=hashlib.sha256(data).hexdigest(),
    data_rows=len(records)-1,columns=width,orientation=orientation,coverage=coverage,
    no_effect_analysis=True,cell_coordinates_not_inherited_from_compartment_average=True))
  print(json.dumps(dict(file=name,download_verified=True,orientation=orientation,required34_measured=len(measured))),flush=True)
 receipt=dict(status='ALL4_SMALL_AUTHOR_PUBLIC_INPUTS_GIT_BLOB_VERIFIED_AND_FIXED_PANEL_COVERAGE_CHECKED',finished_utc=datetime.now(timezone.utc).isoformat(),
  source_study='https://www.nature.com/articles/s41588-025-02351-7',source_repository='https://github.com/tznaung/NSCLC_SpatialOmics',
  source_tree_revision=REVISION,total_bytes=sum(x['bytes'] for x in rows),files=rows,no_patient_rows_in_receipt=True,
  single_cell_contact_not_measured=True,clinical_or_spatial_model_not_started=True)
 (a.output/'public_source_measurement_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);main(p.parse_args())
