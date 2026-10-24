"""Independently cross-check saved large-model 100 K hold outputs, not equilibrium."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
from openmm import XmlSerializer,unit,MonteCarloBarostat,CustomExternalForce
try:
    from .dcd_diagnostics import frames
except ImportError:
    from dcd_diagnostics import frames

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def validate(run):
    run=Path(run);r=json.loads((run/'result.json').read_text());cfg=json.loads((run/'config.json').read_text())
    if r['steps']!=2000 or r['time_ps']!=2 or r['heating_to_310K_completed'] or r['production_ready'] or r['equilibrium_established']:raise ValueError('Incorrect completion/interpretation status')
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    if len(rows)!=81 or [int(x['step']) for x in rows]!=list(range(0,2001,25)):raise ValueError('Observation gaps or step mismatch')
    for x in rows:
        step=int(x['step']);target=100
        if int(x['target_K'])!=target or abs(float(x['time_ps'])-step*.001)>1e-12:raise ValueError('Schedule or time mismatch')
        if not all(math.isfinite(float(v)) for v in x.values()):raise ValueError('Nonfinite record')
    if max(float(x['volume_nm3']) for x in rows)-min(float(x['volume_nm3']) for x in rows)>1e-8:raise ValueError('Fixed-box volume changed')
    system=XmlSerializer.deserialize((run/'system.xml').read_text());state=XmlSerializer.deserialize((run/'state.xml').read_text());integrator=XmlSerializer.deserialize((run/'integrator-final.xml').read_text())
    if system.getNumParticles()!=577890 or any(isinstance(f,MonteCarloBarostat) for f in system.getForces()):raise ValueError('Wrong particle count or barostat')
    restraint=next(f for f in system.getForces() if isinstance(f,CustomExternalForce))
    if restraint.getNumParticles()!=31974 or any(abs(restraint.getParticleParameters(i)[1][0]-1000)>1e-9 for i in range(restraint.getNumParticles())):raise ValueError('Restraint support changed')
    if abs(integrator.getStepSize().value_in_unit(unit.femtosecond)-1)>1e-12 or integrator.getTemperature().value_in_unit(unit.kelvin)!=100 or abs(integrator.getFriction().value_in_unit(unit.picosecond**-1)-5)>1e-12:raise ValueError('Final integrator differs')
    if abs(state.getTime().value_in_unit(unit.picosecond)-3)>1e-10 or state.getStepCount()!=3000:raise ValueError('Saved state time/step differs')
    xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer);vel=state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)
    if not np.isfinite(xyz).all() or not np.isfinite(vel).all():raise ValueError('Nonfinite final coordinates/velocities')
    box=np.array([v.value_in_unit(unit.nanometer) for v in state.getPeriodicBoxVectors()]);volume=float(np.linalg.det(box))
    if abs(volume-float(rows[-1]['volume_nm3']))>1e-6:raise ValueError('CSV and restart cell differ')
    count=0
    for step,p,cell in frames(run/'trajectory.dcd'):count+=1
    if count!=20 or step!=2000 or step+cfg['starting_context_step']!=state.getStepCount() or len(p)!=577890 or not np.allclose(cell,np.diag(box),atol=1e-7):raise ValueError('DCD count, atom identity or cell mismatch')
    delta=p-xyz;delta-=cell*np.round(delta/cell);mismatch=float(np.linalg.norm(delta,axis=1).max())
    if mismatch>1e-5:raise ValueError('DCD and restart coordinates disagree')
    source=resolve_project_path(cfg['source_run']);same_system=hashlib.sha256((run/'system.xml').read_bytes()).hexdigest()==hashlib.sha256((source/'system.xml').read_bytes()).hexdigest()
    if not same_system:raise ValueError('Force field/restraint XML changed during hold')
    initial=XmlSerializer.deserialize((run/'initial-state.xml').read_text());previous=XmlSerializer.deserialize((source/'state.xml').read_text())
    for method,units in [('getPositions',unit.nanometer),('getVelocities',unit.nanometer/unit.picosecond)]:
        a=getattr(initial,method)(asNumpy=True).value_in_unit(units);b=getattr(previous,method)(asNumpy=True).value_in_unit(units)
        if not np.allclose(a,b,atol=1e-10,rtol=0):raise ValueError('Hold initial state differs from previous final state')
    if not cfg['checkpoint_XML_coordinates_and_velocities_match']:raise ValueError('Checkpoint continuity was not checked')
    check={'status':'independent_100K_hold_output_checks_passed','atoms':577890,'steps':2000,'time_ps':2,'cumulative_context_steps':3000,'cumulative_context_time_ps':3,'DCD_frames':count,'DCD_last_segment_local_step':step,'DCD_context_step_offset':cfg['starting_context_step'],'DCD_header_time_is_segment_local':True,'observations':len(rows),'barostats':0,'restrained_heavy_atoms':31974,'source_and_hold_system_XML_identical':same_system,'initial_coordinates_and_velocities_match_previous_final':True,'maximum_DCD_restart_minimum_image_coordinate_difference_nm':mismatch,'saved_velocities_finite':True,'final_integrator_target_K':100,'temperature_equilibrium_established':False,'interface_equilibrium_established':False,'production_ready':False}
    (run/'validation.json').write_text(json.dumps(check,indent=2));print(json.dumps(check,indent=2));return check

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args();validate(a.run)
