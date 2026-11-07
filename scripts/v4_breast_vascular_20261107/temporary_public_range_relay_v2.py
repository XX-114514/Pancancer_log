"""Finite loopback-only relay for two explicitly allowed public research archives.
SSH identity/host and remote Python command are supplied externally. The relay
ends with that finite command; it does not modify proxy settings or expose a
listener on an external interface. All HTTP requests must be exact byte ranges.
"""
import os,subprocess,threading,urllib.request,http.server,re,json
from pathlib import Path
ALLOWED={'/anndata_objects.zip':'https://zenodo.org/records/18917100/files/anndata_objects.zip?download=1','/cell_type_mapping.zip':'https://zenodo.org/records/18917100/files/cell_type_mapping.zip?download=1'}
class Handler(http.server.BaseHTTPRequestHandler):
 def do_GET(self):
  if self.path not in ALLOWED:self.send_error(404);return
  ran=self.headers.get('Range','')
  if not re.fullmatch(r'bytes=\d+-\d+',ran):self.send_error(400);return
  start,end=map(int,ran[6:].split('-'))
  if end<start or end-start>4*1024**2:self.send_error(400);return
  try:
   req=urllib.request.Request(ALLOWED[self.path],headers={'Range':ran,'User-Agent':'PanCancer-public-research/1.0'})
   with urllib.request.urlopen(req,timeout=90) as r:
    if r.status!=206 or not r.headers['Content-Range'].startswith(f'bytes {start}-{end}/'):raise RuntimeError('Range not honored')
    b=r.read(end-start+2);assert len(b)==end-start+1
    self.send_response(206);self.send_header('Content-Range',r.headers['Content-Range']);self.send_header('Content-Length',str(len(b)));self.end_headers();self.wfile.write(b)
  except (BrokenPipeError,ConnectionResetError):pass
  except Exception:self.send_error(502)
 def log_message(self,*args):pass
if __name__=='__main__':
 server=http.server.ThreadingHTTPServer(('127.0.0.1',0),Handler);threading.Thread(target=server.serve_forever,daemon=True).start()
 local_port=server.server_address[1];remote_port=int(os.environ['RELAY_REMOTE_PORT']);host=os.environ['RESEARCH_SSH_HOST'];command=os.environ['FINITE_REMOTE_COMMAND'];w=Path(os.environ['RELAY_WORK']);assert command
 argv=['ssh','-T','-C','-o','ControlPath=none','-o','ControlMaster=no','-o','BatchMode=yes','-o','UpdateHostKeys=no','-o','ExitOnForwardFailure=yes','-o','ConnectTimeout=20','-o','ServerAliveInterval=20','-o','ServerAliveCountMax=3','-R',f'127.0.0.1:{remote_port}:127.0.0.1:{local_port}',host,command]
 (w/(os.environ.get('RELAY_LABEL','relay')+'_transport_PRIVATE.json')).write_text(json.dumps({'listener':'loopback_only','allowed_public_archives':list(ALLOWED),'scope':'one finite remote acquisition command; no persistent proxy or SSH changes'},indent=2))
 try:
  with (w/(os.environ.get('RELAY_LABEL','relay')+'_stdout_PRIVATE.log')).open('ab') as out,(w/(os.environ.get('RELAY_LABEL','relay')+'_stderr_PRIVATE.log')).open('ab') as err:
   p=subprocess.run(argv,stdout=out,stderr=err,timeout=18*3600)
  (w/(os.environ.get('RELAY_LABEL','relay')+'_terminal_PRIVATE.json')).write_text(json.dumps({'exit_code':p.returncode,'finite_command_ended':True},indent=2))
 finally:server.shutdown();server.server_close()
