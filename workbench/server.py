"""Local, same-origin workbench API; compute workers are separate processes."""
import argparse
import json
import mimetypes
import os
from pathlib import Path
import secrets
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.parse import urlsplit,unquote
from .config import ROOT,MODELS,availability,make_config
from .manager import Manager,BusyError


def application(manager):
 token=secrets.token_urlsafe(32)
 class Handler(BaseHTTPRequestHandler):
  def send(self,status,data,ctype='application/json; charset=utf-8'):
   body=json.dumps(data,ensure_ascii=False).encode() if isinstance(data,(dict,list)) else data
   self.send_response(status);self.send_header('Content-Type',ctype);self.send_header('Content-Length',str(len(body)));self.send_header('X-Content-Type-Options','nosniff');self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(body)
  def do_GET(self):
   try:
    parts=unquote(urlsplit(self.path).path).strip('/').split('/')
    if parts==['api','capabilities']:
     return self.send(200,{'token':token,'runtime':manager.runtime,'models':{k:{**v,**availability(k,manager.root)} for k,v in MODELS.items()},'unavailable':['PTFE','PVDF α','PVDF β','PVDF γ']})
    if parts==['api','jobs']:return self.send(200,manager.list())
    if len(parts)>=3 and parts[:2]==['api','jobs']:
     job=manager.directory(parts[2])
     if len(parts)==3:return self.send(200,manager.get(parts[2]))
     if parts[3:]==['observations']:
      rows=[];p=job/'observations.jsonl'
      if p.exists():
       for line in p.read_text().splitlines():
        try:rows.append(json.loads(line))
        except ValueError:break
      return self.send(200,rows)
     if parts[3]=='files' and len(parts)>=5:
      p=(job/('/'.join(parts[4:]))).resolve()
      if not p.is_relative_to(job.resolve()) or p.is_symlink() or any(x.startswith('.') for x in parts[4:]):raise ValueError('Invalid artifact path')
      return self.file(p)
     raise FileNotFoundError()
    if parts and parts[0]=='api':raise FileNotFoundError()
    # Serve public UI/code/docs/artifacts only; raw inputs and hidden files stay private.
    name='/'.join(parts) or 'index.html';p=(manager.root/name).resolve()
    if not p.is_relative_to(manager.root) or any(x.startswith('.') for x in parts) or parts[0] in {'runs','data','datasets','workbench','md','scripts','tests','node_modules'}:raise FileNotFoundError()
    return self.file(p)
   except FileNotFoundError:self.send(404,{'error':'Not found'})
   except (ValueError,OSError) as e:self.send(400,{'error':str(e)})
  def file(self,p):
   if not p.is_file():raise FileNotFoundError()
   kind=mimetypes.guess_type(p.name)[0] or 'application/octet-stream'
   self.send_response(200);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(p.stat().st_size));self.send_header('X-Content-Type-Options','nosniff');self.end_headers()
   with p.open('rb') as f:
    while block:=f.read(1024**2):self.wfile.write(block)
  def do_POST(self):
   try:
    if self.headers.get('X-Workbench-Token')!=token:return self.send(403,{'error':'Missing workbench token'})
    origin=self.headers.get('Origin')
    if origin and urlsplit(origin).netloc!=self.headers.get('Host'):return self.send(403,{'error':'Same-origin requests required'})
    if self.headers.get_content_type()!='application/json':return self.send(415,{'error':'JSON required'})
    n=int(self.headers.get('Content-Length','0'))
    if not 0<n<=8192:return self.send(400,{'error':'Invalid request size'})
    req=json.loads(self.rfile.read(n));parts=urlsplit(self.path).path.strip('/').split('/')
    if parts==['api','preview']:return self.send(200,make_config(req))
    if parts==['api','jobs']:return self.send(201,manager.start(req))
    if len(parts)==4 and parts[:2]==['api','jobs'] and parts[3] in {'pause','resume','cancel'}:return self.send(200,manager.control(parts[2],parts[3]))
    self.send(404,{'error':'Not found'})
   except BusyError as e:self.send(409,{'error':str(e)})
   except FileNotFoundError:self.send(404,{'error':'Not found'})
   except (ValueError,OSError) as e:self.send(400,{'error':str(e)})
 return Handler


def main():
 p=argparse.ArgumentParser();p.add_argument('--host',default='127.0.0.1');p.add_argument('--port',type=int,default=int(os.environ.get('PORT','8000')));a=p.parse_args()
 server=ThreadingHTTPServer((a.host,a.port),application(Manager()));print(f'Workbench: http://{a.host}:{server.server_port}',flush=True)
 try:server.serve_forever()
 except KeyboardInterrupt:pass
 finally:server.server_close()

if __name__=='__main__':main()
