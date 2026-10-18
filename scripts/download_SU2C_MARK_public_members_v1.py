#!/usr/bin/env python3
"""Acquire selected open author ZIP members with HTTP ranges and per-member CRC/SHA."""
import argparse,hashlib,io,json,urllib.request,zipfile
from pathlib import Path
from datetime import datetime,timezone

SOURCE='https://zenodo.org/records/11179623/files/Source%20Data%203.zip?download=1'
MEMBERS=['Source Data/Clinical/Table_S1_Clinical_Annotations.xlsx',
 'Source Data/Clinical/SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt',
 'Source Data/RNA/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct']
class RemoteRange(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.pos=0;self.network_bytes=0;self.cache={};self.block=4*1024**2
  with self.request(0,1) as r:
   assert r.status==206,'Public server must support range requests'
   self.size=int(r.headers['Content-Range'].split('/')[-1]);r.read();self.network_bytes+=1
 def request(self,start,n):
  return urllib.request.urlopen(urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{start+n-1}','User-Agent':'PanCancer-public-research/1.0'}),timeout=60)
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos);parts=[]
  while n>0:
   idx=self.pos//self.block;offset=self.pos%self.block
   if idx not in self.cache:
    start=idx*self.block;size=min(self.block,self.size-start)
    with self.request(start,size) as r:
     assert r.status==206;data=r.read(size+1);assert len(data)==size
    if len(self.cache)>=4:self.cache.pop(next(iter(self.cache)))
    self.cache[idx]=data;self.network_bytes+=len(data);assert self.network_bytes<100*1024**2
   data=self.cache[idx];take=min(n,len(data)-offset);parts.append(data[offset:offset+take]);self.pos+=take;n-=take
  return b''.join(parts)

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700,parents=True)
 reader=RemoteRange(SOURCE);records=[]
 with zipfile.ZipFile(reader) as z:
  total=sum(z.getinfo(n).file_size for n in MEMBERS);assert total<50*1024**2
  for name in MEMBERS:
   info=z.getinfo(name);out=a.output/Path(name).name;h=hashlib.sha256();size=0
   with z.open(name) as src,out.open('xb') as dst:
    for chunk in iter(lambda:src.read(1024**2),b''):
     dst.write(chunk);h.update(chunk);size+=len(chunk)
   assert size==info.file_size
   records.append(dict(name=out.name,archive_member=name,bytes=size,sha256=h.hexdigest(),source_CRC32=info.CRC,full_member_CRC_verified_by_zipfile=True))
   print(json.dumps(dict(phase='public_member_downloaded',name=out.name,bytes=size)),flush=True)
 result=dict(status='SELECTED_OPEN_AUTHOR_MEMBERS_DOWNLOADED_CRC_AND_SHA_VERIFIED',finished_utc=datetime.now(timezone.utc).isoformat(),source_record='https://zenodo.org/records/11179623',source_url=SOURCE,source_archive_md5_published='9bf9af4fb9ae1191ea9562a20f2c1f56',whole_archive_MD5_not_verified=True,archive_bytes=reader.size,network_bytes=reader.network_bytes,members=records,controlled_raw_data_not_requested=True,no_outcome_analysis=True)
 (a.output/'download_receipt_public.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
