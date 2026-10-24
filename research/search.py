import urllib.request, urllib.parse, json
from pathlib import Path
out=Path(__file__).resolve().parent
queries=['"carboxybetaine" AND "Jiang" AND FIRST_PDATE:[2016-01-01 TO 2026-12-31]','"CBAA" AND "Jiang"','"factor XII" AND structure AND zymogen']
for i,q in enumerate(queries):
 u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','pageSize':30,'resultType':'core'})
 d=json.load(urllib.request.urlopen(u,timeout=30)); (out/f'search-{i}.json').write_text(json.dumps(d,indent=2))
 print('\nQUERY',q,'HITS',d.get('hitCount'))
 for x in d.get('resultList',{}).get('result',[]):
  print(json.dumps({k:x.get(k) for k in ['id','pmcid','title','authorString','doi','pubYear','abstractText']},ensure_ascii=False))
