"""Immutable baseline + bounded incremental ZIPs for a complete indexed snapshot."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import uuid
import zipfile

FORMAT='codex-results-layers-v1'
ARCHIVE_NAME=re.compile(r'(results-complete|results-delta-[0-9]{4}-[0-9a-f]{12})\.zip')
MAX_LAYER_BYTES=1_250_000_000


def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while chunk:=f.read(8*1024**2):h.update(chunk)
 return h.hexdigest()


def safe_path(value):
 if not isinstance(value,str):raise ValueError('Unsafe snapshot path')
 p=PurePosixPath(value)
 if p.is_absolute() or '..' in p.parts or '\\' in value or not p.parts or p.as_posix()!=value:raise ValueError('Unsafe snapshot path')
 return p


def atomic(path,data):
 path=Path(path);temp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.part');temp.write_text(json.dumps(data,ensure_ascii=False,indent=2));temp.replace(path)


def load_index(path):
 return validate_index(json.loads(Path(path).read_text()))


def validate_index(index):
 if index.get('format')!=FORMAT:raise ValueError('Unknown archive-index format')
 names=set()
 for a in index['archives']:
  if not ARCHIVE_NAME.fullmatch(a['name']) or a['name'] in names or not re.fullmatch('[0-9a-f]{64}',a['sha256']) or type(a['bytes'])is not int or a['bytes']<0:raise ValueError('Invalid archive identity')
  names.add(a['name'])
 seen=set()
 for f in index['files']:
  safe_path(f['path'])
  if f['path'] in seen or f['archive'] not in names or type(f['bytes'])is not int or f['bytes']<0 or not re.fullmatch('[0-9a-f]{64}',f['sha256']):raise ValueError('Invalid snapshot file entry')
  seen.add(f['path'])
 return index


def initialize(datasets,expected_sha256=None):
 datasets=Path(datasets);path=datasets/'results-index.json'
 if path.exists():return load_index(path)
 baseline=datasets/'results-complete.zip';sha=digest(baseline)
 if expected_sha256 and sha!=expected_sha256:raise ValueError('Baseline archive checksum differs')
 with zipfile.ZipFile(baseline) as z:files=json.loads(z.read('FILES.json'))['files']
 index={'format':FORMAT,'generated_utc':datetime.now(timezone.utc).isoformat(),'archives':[{'name':baseline.name,'bytes':baseline.stat().st_size,'sha256':sha}],'files':[{**f,'archive':baseline.name} for f in files],'note':'Complete latest snapshot selected by file SHA-256. Immutable baseline and delta archives also preserve previous versions. Restore data/runs through this index, not by overwriting files in archive order.'}
 validate_index(index);atomic(path,index);return index


def publish_snapshot(full_snapshot,datasets):
 full_snapshot=Path(full_snapshot);datasets=Path(datasets);index=initialize(datasets);previous={f['path']:f for f in index['files']};archives=list(index['archives']);new=[]
 with zipfile.ZipFile(full_snapshot) as source:
  manifest=json.loads(source.read('FILES.json'));files=manifest['files'];seen=set()
  for f in files:
   safe_path(f['path'])
   if f['path'] in seen:raise ValueError('Duplicate snapshot member')
   seen.add(f['path'])
   if f['path'] not in previous or any(previous[f['path']][k]!=f[k] for k in ['bytes','sha256']):new.append(f)
  batches=[];batch=[];size=0
  for f in new:
   if f['bytes']>MAX_LAYER_BYTES:raise ValueError('Individual snapshot member exceeds bounded ZIP budget')
   if batch and size+f['bytes']>MAX_LAYER_BYTES:batches.append(batch);batch=[];size=0
   batch.append(f);size+=f['bytes']
  if batch:batches.append(batch)
  chosen={f['path']:previous[f['path']]['archive'] for f in files if f['path'] in previous and all(previous[f['path']][k]==f[k] for k in ['bytes','sha256'])}
  for batch in batches:
   temp=datasets/('delta-'+uuid.uuid4().hex+'.zip.part')
   try:
    with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as target:
     for entry in batch:
      h=hashlib.sha256();length=0
      with source.open(entry['path']) as inp,target.open(entry['path'],'w',force_zip64=True) as out:
       while chunk:=inp.read(8*1024**2):length+=len(chunk);h.update(chunk);out.write(chunk)
      if length!=entry['bytes'] or h.hexdigest()!=entry['sha256']:raise ValueError('Delta member integrity mismatch')
     target.writestr('FILES.json',json.dumps({'files':batch,'generated_utc':manifest.get('generated_utc'),'scope':'Incremental layer; use datasets/results-index.json for a complete snapshot.'},ensure_ascii=False,indent=2))
    sha=digest(temp);name=f'results-delta-{len(archives):04d}-{sha[:12]}.zip';dest=datasets/name
    if dest.exists():raise ValueError('Refusing to overwrite an immutable delta')
    temp.replace(dest);archives.append({'name':name,'bytes':dest.stat().st_size,'sha256':sha})
    for f in batch:chosen[f['path']]=name
   finally:temp.unlink(missing_ok=True)
 updated={**index,'generated_utc':manifest.get('generated_utc',datetime.now(timezone.utc).isoformat()),'archives':archives,'files':[{**f,'archive':chosen[f['path']]} for f in files]}
 validate_index(updated);atomic(datasets/'results-index.json',updated);return updated
