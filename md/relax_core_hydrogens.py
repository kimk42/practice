"""Pre-relax generated hydrogens/caps with observed heavy atoms strongly restrained.
Unsolvated, truncated reaction-field setup is for local preparation only, NOT adsorption MD.
"""
import argparse,json,math,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import unit,Context,VerletIntegrator,Platform,LocalEnergyMinimizer,MinimizationReporter,CustomExternalForce
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,PDBFile,CutoffNonPeriodic,element,HBonds


def close_contacts(topology,xyz,cutoff=.12):
    atoms=list(topology.atoms());bonds={frozenset((a.index,b.index)) for a,b in topology.bonds()}
    pairs=cKDTree(xyz).query_pairs(cutoff)
    return [(i,j) for i,j in pairs if frozenset((i,j)) not in bonds]

def relax(bundle,iterations=100):
    root=Path(bundle);out=root/'local-relaxation-restrained';out.mkdir(exist_ok=False)
    manifest=json.loads((root/'manifest.json').read_text())
    p=AmberPrmtopFile(str(root/'core.prmtop'));c=AmberInpcrdFile(str(root/'core.inpcrd'));original=np.array(c.positions.value_in_unit(unit.nanometer))
    system=p.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1*unit.nanometer,constraints=HBonds,rigidWater=False)
    fixed=[]
    restraint=CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addGlobalParameter('k',100000)
    for name in ('x0','y0','z0'):restraint.addPerParticleParameter(name)
    for a in p.topology.atoms():
        # Source heavy atoms and calcium fixed; artificial caps and hydrogens may relax.
        if a.element!=element.hydrogen and not manifest['residue_mapping'][a.residue.index]['artificial_cap']:
            restraint.addParticle(a.index,original[a.index]);fixed.append(a.index)
    system.addForce(restraint)
    integrator=VerletIntegrator(.001*unit.picoseconds)
    context=Context(system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});context.setPositions(c.positions)
    initial=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    class Progress(MinimizationReporter):
        def report(self,iteration,x,grad,args):
            if iteration%20==0:print('Local H/cap relaxation',iteration,args['system energy'],flush=True)
            return False
    LocalEnergyMinimizer.minimize(context,10*unit.kilojoule_per_mole/unit.nanometer,iterations,Progress())
    state=context.getState(getEnergy=True,getPositions=True);xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy) or not np.isfinite(xyz).all():raise ValueError('Non-finite local relaxation')
    displacement=float(np.max(np.linalg.norm(xyz[fixed]-original[fixed],axis=1)))
    if displacement>.1:raise ValueError('Observed atom displacement exceeds preparation bound')
    with (out/'core-relaxed.pdb').open('w') as f:PDBFile.writeFile(p.topology,state.getPositions(),f)
    # Amber coordinates retain exact topology atom order; no new topology is inferred from PDB.
    flat=(xyz*10).ravel()
    with (out/'core-relaxed.inpcrd').open('w') as f:
        f.write('3GHG capped observed core: local H/cap preparation only\n');f.write(f'{len(xyz):6d}\n')
        for i in range(0,len(flat),6):f.write(''.join(f'{v:12.7f}' for v in flat[i:i+6])+'\n')
    after=close_contacts(p.topology,xyz);before=close_contacts(p.topology,original)
    cap_charges={}
    from openmm import NonbondedForce
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    for r in p.topology.residues():
        if manifest['residue_mapping'][r.index]['artificial_cap']:
            cap_charges[str(r.index+1)]=sum(nb.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge) for a in r.atoms())
    if max(abs(q) for q in cap_charges.values())>1e-5:raise ValueError('Artificial cap not neutral')
    report={'status':'local_preparation_completed_not_equilibrated','initial_diagnostic_energy_kJ_mol':initial,'final_diagnostic_energy_kJ_mol':energy,
            'restrained_observed_heavy_atoms':len(fixed),'heavy_restraint_k_kJ_mol_nm2':100000,'max_restrained_atom_displacement_nm':displacement,'free_particles':'generated H and artificial caps; observed heavy atoms strongly restrained, X-H bonds constrained','cap_charges_e':cap_charges,
            'nonbonded_pairs_under_1_2A_before':len(before),'nonbonded_pairs_under_1_2A_after':len(after),
            'iterations_limit':iterations,'convergence':'not assessed; iteration-limited preparation','setup':'Unsolvated ff14SB with heavy positional restraints and HBonds constraints, nonperiodic 1-nm reaction-field cutoff; diagnostic H/cap preparation only, not a solvent model for MD',
            'source_prmtop_sha256':hashlib.sha256((root/'core.prmtop').read_bytes()).hexdigest(),'relaxed_inpcrd_sha256':hashlib.sha256((out/'core-relaxed.inpcrd').read_bytes()).hexdigest(),
            'interpretation':'No physical binding energy, adsorption, equilibrium, or material ranking is inferred.'}
    (out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle');p.add_argument('--iterations',type=int,default=100);a=p.parse_args();relax(a.bundle,a.iterations)
