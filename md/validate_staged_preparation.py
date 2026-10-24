"""Check execution integrity, ensembles and saved periodic cells; not an equilibrium test."""
import argparse,csv,json,math,struct,hashlib
from pathlib import Path
import numpy as np
from openmm import XmlSerializer,unit,MonteCarloBarostat,CustomExternalForce


def validate(run):
    root=Path(run);config=json.loads((root/'config.json').read_text());result=json.loads((root/'result.json').read_text())
    with (root/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    expected=sum(n for stage in config['stages'] for _,n in stage['schedule'])
    if result['steps']!=expected or int(rows[-1]['cumulative_step'])!=expected:raise ValueError('Incomplete executed steps')
    if abs(result['time_ps']-expected*.002)>1e-10:raise ValueError('Time mismatch')
    if result['equilibration_established'] or result['production_ready']:raise ValueError('Short pilot incorrectly labeled ready')
    for row in rows:
        if not all(math.isfinite(float(v)) for k,v in row.items() if k!='stage'):raise ValueError('Non-finite observation')
    checks=[]
    for spec,stage_result in zip(config['stages'],result['stages']):
        directory=root/spec['name'];system=XmlSerializer.deserialize((directory/'system.xml').read_text());state=XmlSerializer.deserialize((directory/'state.xml').read_text())
        coords=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        if system.getNumParticles()!=43288 or len(coords)!=43288:raise ValueError('Particle count mismatch')
        bars=sum(isinstance(f,MonteCarloBarostat) for f in system.getForces());restraints=sum(f.getNumParticles() for f in system.getForces() if isinstance(f,CustomExternalForce))
        if bars!=int(spec['ensemble']=='NPT'):raise ValueError('Barostat/ensemble mismatch')
        if spec['ensemble']=='NPT' and restraints!=0:raise ValueError('NPT unexpectedly retains Cartesian restraints')
        subset=[r for r in rows if r['stage']==spec['name']];vol=[float(r['volume_nm3']) for r in subset]
        if spec['ensemble']=='NVT' and max(vol)-min(vol)>1e-7:raise ValueError('NVT volume changed')
        if spec['ensemble']=='NPT' and max(vol)-min(vol)<1e-6:raise ValueError('No barostat volume response observed')
        if stage_result['executed_steps']!=sum(n for _,n in spec['schedule']):raise ValueError('Stage incomplete')
        error=0
        for i in range(system.getNumConstraints()):
            a,b,d=system.getConstraintParameters(i);target=d.value_in_unit(unit.nanometer)
            error=max(error,abs(np.linalg.norm(coords[a]-coords[b])-target)/target)
        if error>.001:raise ValueError('Constraint error too large')
        with (directory/'trajectory.dcd').open('rb') as f:header=f.read(12)
        if header[4:8]!=b'CORD':raise ValueError('Invalid DCD header')
        frames=struct.unpack('<i',header[8:12])[0]
        if frames!=stage_result['executed_steps']//100:raise ValueError('DCD frame count mismatch')
        vectors=np.array([v.value_in_unit(unit.nanometer) for v in state.getPeriodicBoxVectors()]);volume=float(np.linalg.det(vectors))
        if abs(volume-vol[-1])>1e-6:raise ValueError('State box differs from recorded volume')
        # For the initial implementation the active process predates the topology-box writer fix.
        # Rebuild CRYST1 from the authoritative saved State, preserving all atom coordinates.
        lengths=np.linalg.norm(vectors,axis=1)*10
        angles=[np.degrees(np.arccos(np.clip(np.dot(vectors[a],vectors[b])/(np.linalg.norm(vectors[a])*np.linalg.norm(vectors[b])),-1,1))) for a,b in [(1,2),(0,2),(0,1)]]
        crystal=f'CRYST1{lengths[0]:9.3f}{lengths[1]:9.3f}{lengths[2]:9.3f}{angles[0]:7.2f}{angles[1]:7.2f}{angles[2]:7.2f} P 1           1 '
        pdb=directory/'final.pdb';lines=pdb.read_text().splitlines();old=next((l for l in lines if l.startswith('CRYST1')),None)
        lines=[crystal if l.startswith('CRYST1') else l for l in lines]
        if old is None:lines.insert(1,crystal)
        pdb.write_text('\n'.join(lines)+'\n')
        checks.append({'stage':spec['name'],'particles':43288,'barostats':bars,'restrained_particles':restraints,'max_relative_constraint_error':float(error),'DCD_frames':frames,'PDB_CRYST1_from_saved_state':True,'PDB_cell_header_updated':old!=crystal,'state_sha256':hashlib.sha256((directory/'state.xml').read_bytes()).hexdigest()})
    report={'status':'execution_integrity_checks_passed','scientific_equilibrium_established':False,'stages':checks,'total_steps_verified':expected,'recorded_samples':len(rows)}
    (root/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run');a=p.parse_args();validate(a.run)
