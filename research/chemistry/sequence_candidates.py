"""Illustrative exact-count shuffled sequences, not measured chain architecture."""
import argparse, collections, json, random
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--dp',type=int,default=100);p.add_argument('--replicates',type=int,default=5);a=p.parse_args()
if a.dp<=0 or a.dp%10:raise SystemExit('DP must be a positive multiple of 10 for exact 7:2:1 counts')
if a.replicates<=0:raise SystemExit('replicates must be positive')
counts={'CBAA1':7*a.dp//10,'BMA':2*a.dp//10,'BPAA':a.dp//10}
chains=[]
for seed in range(a.replicates):
 seq=[m for m,n in counts.items() for _ in range(n)];random.Random(seed).shuffle(seq)
 assert dict(collections.Counter(seq))==counts
 chains.append({'seed':seed,'sequence':seq,'bpaa_indices_1based':[i+1 for i,x in enumerate(seq) if x=='BPAA']})
d={'status':'sequence candidates only; no atoms, force field or crosslinks',
 'assumptions':['7:2:1 interpreted as incorporated molar composition for this example only; experimental basis pending',
 'DP is illustrative and not measured', 'uniform shuffle conditional on exact counts is a modeling choice, not polymerization kinetics'],
 'dp':a.dp,'counts':counts,'chains':chains}
f=Path(__file__).with_name(f'sequences-dp{a.dp}.json');f.write_text(json.dumps(d,indent=2)+'\n');print('Wrote',f,'with counts',counts,'for',len(chains),'seeds')
