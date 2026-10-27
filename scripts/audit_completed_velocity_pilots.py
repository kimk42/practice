"""Read-only verification of18 completed short pilots;no simulation restarts."""
import json,hashlib,re
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    index=json.loads((ROOT/'runs/velocity-pilots/result.json').read_text());records=index['records'];seen=set();rows=[]
    if len(records)!=18:raise ValueError('Expected18 registered pilots')
    for row in records:
        key=(row['case'],row['seed'])
        if key in seen:raise ValueError('Duplicate case/seed')
        seen.add(key);result_path=ROOT/row['source_result'];out=result_path.parent;result=json.loads(result_path.read_text())
        if result.get('duration_ps')!=.2:raise ValueError('Duration mismatch')
        mdp=(out/'md.mdp').read_text();seed=int(re.search(r'gen-seed\s*=\s*(\d+)',mdp).group(1))
        if seed!=row['seed']:raise ValueError('Seed mismatch')
        log=(out/'md.log').read_text()
        if 'Finished mdrun' not in log or re.search('LINCS WARNING|constraint error',log,re.I):raise ValueError('Incomplete/constraint-failed pilot')
        data=np.array([[float(x) for x in line.split()] for line in (out/'energies.xvg').read_text().splitlines() if line and not line.startswith(('#','@'))])
        if not np.isfinite(data).all() or data[-1,0]-data[0,0]<.199:raise ValueError('Energy duration/finite gate failed')
        rows.append({'case':row['case'],'seed':seed,'duration_ps':float(data[-1,0]-data[0,0]),'energy_frames':len(data),'endpoint_sha256':hashlib.sha256((out/'md.gro').read_bytes()).hexdigest(),'state':'short numerical checks passed'})
    report={'state':'18 completed numerical pilots independently rechecked','records':rows,'scope':'different velocities/orientations/gap subset,not independent morphology or converged adsorption replicas'}
    (ROOT/'models/validation/velocity-pilots-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(report['state'])
if __name__=='__main__':main()
