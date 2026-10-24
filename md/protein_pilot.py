"""Actual uploaded fibrinogen gamma fragment: Amber14/GBn2 preparation and short CPU pilot.
This is a protein-only preparation check, not a polymer adsorption calculation.
"""
import argparse,csv,hashlib,json,math,os,time
from pathlib import Path
import numpy as np
from openmm import LangevinMiddleIntegrator,Platform,XmlSerializer,MinimizationReporter,unit
from openmm.app import PDBFile,PDBxFile,Modeller,ForceField,Simulation,NoCutoff,HBonds,DCDReporter
from pdbfixer import PDBFixer


def aligned_rmsd(current,reference):
    a=current-current.mean(axis=0);b=reference-reference.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
    return float(np.sqrt(np.mean(np.sum((a@u@correction@vt-b)**2,axis=1))))


def run(source,output,steps=1000):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    source=Path(source)
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    cif=MMCIF2Dict(str(source))
    if cif['_entry.id'][0]!='3GHG':raise ValueError('This preparation is specifically reviewed for uploaded 3GHG')
    pdb=PDBxFile(str(source));mod=Modeller(pdb.topology,pdb.positions)
    # Author numbering checked against label sequence below, rather than assuming it.
    ids={cif['_atom_site.auth_seq_id'][i] for i,ch in enumerate(cif['_atom_site.label_asym_id']) if ch=='C' and cif['_atom_site.label_seq_id'][i] not in ('.','?') and 143<=int(cif['_atom_site.label_seq_id'][i])<=394}
    selected=[r for r in mod.topology.residues() if r.chain.id=='C' and r.id in ids and r.name not in ('HOH','CA')]
    if len(selected)!=252:raise ValueError(f'Expected 252 consecutive residues; got {len(selected)}')
    allowed=set(selected);mod.delete([r for r in mod.topology.residues() if r not in allowed])
    with (out/'observed_fragment.pdb').open('w') as f:PDBFile.writeFile(mod.topology,mod.positions,f)
    fixer=PDBFixer(filename=str(out/'observed_fragment.pdb'))
    fixer.findMissingResidues();fixer.missingResidues={}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:raise ValueError('Nonstandard residues require explicit treatment')
    fixer.findMissingAtoms()
    missing={f'{r.chain.id}:{r.name}{r.id}':[a.name for a in atoms] for r,atoms in fixer.missingAtoms.items()}
    terminals={f'{r.chain.id}:{r.name}{r.id}':names for r,names in fixer.missingTerminals.items()}
    # Do not rebuild absent protein heavy atoms silently; only terminal atoms may be added.
    if missing:raise ValueError(f'Unexpected missing heavy atoms: {missing}')
    fixer.addMissingAtoms(seed=42)
    ff=ForceField('amber14-all.xml','implicit/gbn2.xml')
    mod=Modeller(fixer.topology,fixer.positions)
    mod.addHydrogens(ff,pH=7.4,platform=Platform.getPlatformByName('CPU'))
    system=ff.createSystem(mod.topology,nonbondedMethod=NoCutoff,constraints=None,soluteDielectric=1.0,solventDielectric=78.5)
    (out/'system.xml').write_text(XmlSerializer.serialize(system))
    with (out/'prepared.pdb').open('w') as f:PDBFile.writeFile(mod.topology,mod.positions,f)
    manifest={'status':'prepared','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_entry':'3GHG','chain':'C','label_sequence_range':[143,394],
              'output_residue_mapping':{'note':'PDB export renumbers selected fragment from 1; add 142 to obtain source label sequence index','output_chain':'A','source_chain':'C'},'residues':252,'particles':system.getNumParticles(),'protein_forcefield':'Amber14 (OpenMM amber14-all.xml)','solvent':'GBn2 implicit solvent; solvent dielectric 78.5; no explicit water/ions',
              'pH_for_hydrogen_addition':7.4,'added_terminal_atoms':terminals,'missing_heavy_atoms':missing,
              'disulfides':[f'{a.residue.id}:{a.name}-{b.residue.id}:{b.name}' for a,b in mod.topology.bonds() if a.name==b.name=='SG'],
              'assumptions':['Isolated observed gamma-chain fragment, not intact fibrinogen','Artificial charged peptide termini at fragment boundaries','Other chains, glycans, calcium and knob peptides excluded','No polymer, no adsorption prediction','Hydrogen placement preparation is not guaranteed bitwise reproducible'],
              'minimization_max_iterations':25,'minimization_convergence':'not assessed','steps':steps,'timestep_fs':1,'temperature_K':310,'platform':'CPU','threads':4,'constraints':'none; all bonded terms retained','seed':42}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('Prepared',manifest['particles'],'particles;',manifest['disulfides'],flush=True)
    integrator=LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,1*unit.femtoseconds);integrator.setRandomNumberSeed(42)
    sim=Simulation(mod.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});sim.context.setPositions(mod.positions)
    def state():
        s=sim.context.getState(getEnergy=True,getPositions=True)
        energy=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):raise ValueError('Non-finite energy')
        return s,energy
    _,initial=state();print('Initial potential energy:',initial,flush=True)
    class Progress(MinimizationReporter):
        def report(self,iteration,x,grad,args):
            if iteration%5==0:print('Minimization iteration',iteration,'potential',args['system energy'],flush=True)
            return False
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=25,reporter=Progress())
    minimized,min_energy=state();reference=minimized.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    with (out/'minimized.pdb').open('w') as f:PDBFile.writeFile(mod.topology,minimized.getPositions(),f)
    sim.context.setVelocitiesToTemperature(310*unit.kelvin,42)
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'),100))
    ca=np.array([a.index for a in mod.topology.atoms() if a.name=='CA'])
    dof=3*system.getNumParticles()-system.getNumConstraints()-3
    rows=[];start=time.monotonic()
    with (out/'observables.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['step','time_ps','potential_kJ_mol','temperature_K','CA_RMSD_from_minimized_A'])
        for block in range(0,steps,100):
            sim.step(min(100,steps-block));s,energy=state();xyz=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            kinetic=s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            temp=2*kinetic/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
            row=[sim.currentStep,sim.currentStep*.001,energy,temp,aligned_rmsd(xyz[ca],reference[ca])*10]
            rows.append(row);writer.writerow(row);f.flush();print('Step',sim.currentStep,'energy',round(energy,2),'CA RMSD',round(row[-1],3),'A',flush=True)
    with (out/'final.pdb').open('w') as f:PDBFile.writeFile(mod.topology,s.getPositions(),f)
    result={**manifest,'status':'completed','elapsed_MD_seconds':time.monotonic()-start,'initial_potential_kJ_mol':initial,'minimized_potential_kJ_mol':min_energy,'final_potential_kJ_mol':rows[-1][2],'final_CA_RMSD_from_minimized_A':rows[-1][-1],'interpretation':'Short protein-only preparation pilot. No adsorption, affinity or equilibrium claim.'}
    (out/'result.json').write_text(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output');p.add_argument('--steps',type=int,default=1000)
    a=p.parse_args()
    if not 100<=a.steps<=10000:p.error('steps must be 100..10000')
    print(json.dumps(run(a.source,a.output,a.steps),indent=2))
