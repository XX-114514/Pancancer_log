"""Read official ZIP metadata and selected members using verified HTTP ranges.
Whole-archive MD5 is not certified by a partial download. Every extracted member
is fully read through zipfile (CRC validated), hashed, and bound to archive metadata.
"""
import io,json,urllib.request,zipfile,hashlib,time,argparse
from pathlib import Path
from datetime import datetime,timezone
class HTTPRange(io.RawIOBase):
 def __init__(self,url,budget):
  self.url=url;self.pos=0;self.budget=budget;self.cache={};self.network_bytes=0;self.block=256*1024
  with self.request(0,1) as r:
   self.size=int(r.headers['Content-Range'].split('/')[-1]);assert r.read()==b'P';self.network_bytes+=1
 def request(self,start,n):
  for attempt in range(3):
   try:
    r=urllib.request.urlopen(urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{start+n-1}','User-Agent':'PanCancer-public-research/1.0'}),timeout=300)
    assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {start}-{start+n-1}/'), 'Exact HTTP range not honored'
    return r
   except Exception:
    if attempt==2:raise
    time.sleep(1+attempt)
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):
  self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n
  assert 0<=self.pos<=self.size
  return self.pos
 def read(self,n=-1):
  n=self.size-self.pos if n<0 else min(n,self.size-self.pos);parts=[]
  while n>0:
   idx=self.pos//self.block;offset=self.pos%self.block
   if idx not in self.cache:
    start=idx*self.block;size=min(self.block,self.size-start);assert self.network_bytes+size<=self.budget,'Network budget exceeded'
    for body_attempt in range(3):
     try:
      with self.request(start,size) as r:data=r.read(size+1)
      break
     except Exception:
      if body_attempt==2:raise
      time.sleep(1+body_attempt)
    assert len(data)==size
    if len(self.cache)>=4:self.cache.pop(next(iter(self.cache)))
    self.cache[idx]=data;self.network_bytes+=size
   data=self.cache[idx];take=min(n,len(data)-offset);assert take>0
   parts.append(data[offset:offset+take]);self.pos+=take;n-=take
  return b''.join(parts)
def inventory(url,budget=32*1024**2):
 reader=HTTPRange(url,budget)
 with zipfile.ZipFile(reader) as z:
  members=[{'member':x.filename,'bytes':x.file_size,'compressed_bytes':x.compress_size,'CRC32':x.CRC} for x in z.infolist() if not x.is_dir()]
 return {'archive_bytes':reader.size,'metadata_network_bytes':reader.network_bytes,'members':members}
def download_selected(url,members,out,budget):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);reader=HTTPRange(url,budget);rows=[]
 with zipfile.ZipFile(reader) as z:
  for i,name in enumerate(members):
   info=z.getinfo(name);target=out/f'member_{i:03d}_PRIVATE{Path(name).suffix}';h=hashlib.sha256();size=0
   with z.open(name) as f,target.open('xb') as dst:
    for b in iter(lambda:f.read(1024**2),b''):dst.write(b);h.update(b);size+=len(b)
   assert size==info.file_size
   row={'member':name,'local_name':target.name,'bytes':size,'sha256':h.hexdigest(),'CRC32':info.CRC,'member_CRC_validated':True};rows.append(row)
   (out/f'member_{i:03d}_accepted_PRIVATE.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps({'phase':'member_CRC_SHA_accepted','ordinal':i,'bytes':size}),flush=True)
 receipt={'UTC':datetime.now(timezone.utc).isoformat(),'source_url':url,'archive_bytes':reader.size,'network_bytes':reader.network_bytes,'whole_archive_MD5_verified':False,'members':rows}
 (out/'complete_PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--url',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(inventory(a.url),indent=2)+'\n')
