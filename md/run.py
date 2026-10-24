"""Run an explicitly parameterized protein/polymer system; never invent force fields."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def validate(c):
    surface = c['surface']
    if surface['polymer'] not in ('PTFE', 'PVDF', 'PE'):
        raise ValueError('Unsupported polymer')
    if surface['polymer'] == 'PVDF':
        if surface['phase'] not in ('alpha', 'beta', 'gamma'):
            raise ValueError('PVDF requires alpha/beta/gamma phase')
        if surface['phase'] == 'alpha' and surface['polarization'] != 'nonpolar':
            raise ValueError('This alpha-phase model requires nonpolar packing')
    elif surface['phase'] != 'specified-in-model':
        raise ValueError('Non-PVDF crystal structure must be specified in the model')
    if surface['polarization'] not in ('nonpolar','plus-normal','minus-normal','in-plane','unpoled'):
        raise ValueError('Unsupported polarization')
    if not str(surface['face']).strip():
        raise ValueError('Exposed crystal face is required')
    if not 0 < c['temperature_K'] <= 500 or not 0 < c['timestep_fs'] <= 2:
        raise ValueError('Invalid temperature or timestep')
    if type(c['steps']) is not int or not 1 <= c['steps'] <= 10000000:
        raise ValueError('Invalid step count')
    if type(c['seed']) is not int or not 1 <= c['seed'] < 2147483647:
        raise ValueError('Invalid random seed')
    return c


def run(config_path, bundle, output, platform='CPU'):
    from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, unit
    from openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter
    c = validate(json.loads(Path(config_path).read_text()))
    root = Path(bundle)
    manifest = json.loads((root/'manifest.json').read_text())
    # Require exact model identity: changing UI metadata cannot change physical parameters.
    if manifest['surface'] != c['surface']:
        raise ValueError('Surface configuration differs from the prepared atomic model')
    for key in ('protein_source','polymer_forcefield','protein_forcefield','water_model','preparation_notes'):
        if not manifest.get(key):
            raise ValueError('Missing model provenance: '+key)
    for name in ('system.xml', 'coordinates.pdb'):
        actual = hashlib.sha256((root/name).read_bytes()).hexdigest()
        if manifest['sha256'].get(name) != actual:
            raise ValueError('Model checksum mismatch: '+name)
    pdb = PDBFile(str(root/'coordinates.pdb'))
    system = XmlSerializer.deserialize((root/'system.xml').read_text())
    if system.getNumParticles() != pdb.topology.getNumAtoms():
        raise ValueError('Particle count differs from coordinate topology')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    integrator = LangevinMiddleIntegrator(c['temperature_K']*unit.kelvin, 1/unit.picosecond, c['timestep_fs']*unit.femtoseconds)
    integrator.setRandomNumberSeed(c['seed'])
    kwargs = {'Threads':'2'} if platform == 'CPU' else {}
    sim = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName(platform), kwargs)
    sim.context.setPositions(pdb.positions)
    def energy():
        value = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(value):
            raise ValueError('Non-finite potential energy')
        return value
    initial = energy()
    sim.minimizeEnergy(maxIterations=100)
    sim.context.setVelocitiesToTemperature(c['temperature_K']*unit.kelvin, c['seed'])
    interval = min(100, c['steps'])
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'), interval))
    sim.reporters.append(StateDataReporter(str(out/'energy.csv'), interval, step=True, potentialEnergy=True, temperature=True))
    sim.step(c['steps'])
    final = energy()
    with (out/'final.pdb').open('w') as f:
        PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), f)
    report = {'status':'completed','config':c,'model_manifest':manifest,'platform':platform,'initial_potential_kJ_mol':initial,'final_potential_kJ_mol':final,'interpretation':'Total system potential energy, NOT binding or adsorption free energy.'}
    (out/'result.json').write_text(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('config'); p.add_argument('bundle'); p.add_argument('output')
    p.add_argument('--platform', choices=['CPU','CUDA','OpenCL','Reference'], default='CPU')
    a = p.parse_args()
    print(json.dumps(run(a.config, a.bundle, a.output, a.platform), indent=2))
