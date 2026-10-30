"""Continue each unchanged citrate water candidate from2ps to102ps.

These are descriptive single-solute controls, not the full cPBSzI mixture.
Existing coordinates, velocities, TPR parameters and checkpoints are retained.
"""
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from rdkit import Chem

from pc_length_parameters import ROOT, sha
from run_pc_sqm_checkpointed import write_json

BASE = ROOT/'runs/citrate-solution-102ps-v1'


def table(path):
    return np.array([[float(v) for v in line.split()] for line in path.read_text().splitlines()
                     if line.strip() and not line.startswith(('#', '@'))])


def run(command, out, name, stdin=None, lock_fd=None):
    kwargs = {'pass_fds': (lock_fd,)} if lock_fd is not None else {}
    with (out/(name+'.log')).open('a') as log:
        subprocess.run(['gmx', *command], cwd=out, input=stdin, text=True,
                       stdout=log, stderr=subprocess.STDOUT, check=True, **kwargs)


def analyze(out, source, name, audit, selected_oxygens=None, oxygen_label='carboxylate', maximum_bond_nm=.21):
    target = out/'analysis'
    if (target/'result.json').exists():
        return json.loads((target/'result.json').read_text())
    if target.exists():
        raise RuntimeError('Partial analysis preserved; review before retry')
    logs = list(out.glob('nvt.part*.log'))
    energies = list(out.glob('nvt.part*.edr'))
    trajectories = list(out.glob('nvt.part*.trr'))
    if not (len(logs) == len(energies) == len(trajectories) == 1):
        raise ValueError('Review unexpected continuation output parts')
    text = logs[0].read_text()
    if 'Finished mdrun' not in text or re.search(r'LINCS WARNING|constraint.*warning|Fatal error', text, re.I):
        raise ValueError('Completion or constraint check failed')
    target.mkdir()
    run(['energy', '-f', str(energies[0]), '-o', 'energy.xvg'], target, 'energy', 'Potential\nTemperature\n0\n')
    data = table(target/'energy.xvg')
    if data.shape != (2001, 3) or not np.isfinite(data).all() or not np.allclose(data[[0,-1],0], [2.,102.]):
        raise ValueError('Missing/nonfinite energy frames or incomplete continuation')
    mol = Chem.MolFromMolFile(str(ROOT/f'models/buffer-parameters/{audit["species"]}.mol'), removeHs=False)
    n = mol.GetNumAtoms()
    selected = list(selected_oxygens) if selected_oxygens is not None else []
    if selected_oxygens is None:
        for atom in mol.GetAtoms():
            oxygens = [a for a in atom.GetNeighbors() if a.GetAtomicNum() == 8]
            if atom.GetAtomicNum() == 6 and len(oxygens) == 2 and any(a.GetFormalCharge() == -1 for a in oxygens):
                selected.extend(a.GetIdx()+1 for a in oxygens)
    rows = (source/'solvated.gro').read_text().splitlines()[2:-1]
    groups = {'SOLUTE': list(range(1,n+1)), 'ANION_O': selected,
              'SODIUM': [i+1 for i,s in enumerate(rows) if s[5:10].strip() == 'NAJC'],
              'WATER_O': [i+1 for i,s in enumerate(rows) if s[5:10].strip() == 'SOL' and s[10:15].strip() == 'OW']}
    if not all(groups.values()) or len(groups['SODIUM']) != audit['sodium_ions']:
        raise ValueError('Unexpected chemical/index groups')
    (target/'groups.ndx').write_text(''.join('[ '+k+' ]\n'+' '.join(map(str,v))+'\n' for k,v in groups.items()))
    traj = ['-s', str(out/'nvt.tpr'), '-f', str(trajectories[0]), '-n', 'groups.ndx']
    run(['trjconv', *traj, '-o', 'solute.g96', '-ndec', '9'], target, 'geometry', 'SOLUTE\n')
    lines = (target/'solute.g96').read_text().splitlines()
    xyz = np.array([[[float(v) for v in row.split()] for row in lines[i+1:i+1+n]]
                    for i,line in enumerate(lines) if line == 'POSITIONRED'])
    if xyz.shape != (2001,n,3) or not np.isfinite(xyz).all():
        raise ValueError('Unexpected trajectory coordinates')
    lengths = []
    for bond in mol.GetBonds():
        delta = xyz[:,bond.GetBeginAtomIdx()]-xyz[:,bond.GetEndAtomIdx()]
        delta -= 3*np.rint(delta/3)
        lengths.extend(np.linalg.norm(delta,axis=1).tolist())
    if min(lengths) < .07 or max(lengths) > maximum_bond_nm:
        raise ValueError('Bond geometry outside numerical envelope')
    rdfs = {}
    for window, start, end in [('first',2.,51.95), ('second',52.,102.)]:
        run(['rdf', *traj, '-ref', 'group "ANION_O"', '-sel', 'group "WATER_O"',
             '-b', str(start), '-e', str(end), '-bin', '.01', '-rmax', '1.2',
             '-o', window+'.xvg', '-cn', window+'-cn.xvg'], target, 'rdf-'+window)
        rdf = table(target/(window+'.xvg')); cn = table(target/(window+'-cn.xvg'))
        if not np.isfinite(rdf).all() or not np.isfinite(cn).all():
            raise ValueError('Nonfinite RDF')
        region = rdf[(rdf[:,0]>=.2)&(rdf[:,0]<=.45)]
        peak = region[np.argmax(region[:,1])]
        rdfs[window] = {'saved_time_window_ps':[start,end], 'peak_search_nm':[.2,.45],
                        'peak_r_nm':float(peak[0]),'peak_g_r':float(peak[1]),
                        'water_O_coordination_at_0_35nm':float(np.interp(.35,cn[:,0],cn[:,1]))}
    run(['mindist', *traj, '-od', 'oxygen-na-distance.xvg', '-on', 'oxygen-na-contact.xvg', '-d', '.35'],
        target, 'ion-distance', 'ANION_O\nSODIUM\n')
    distance = table(target/'oxygen-na-distance.xvg')
    if len(distance) != 2001 or not np.isfinite(distance).all():
        raise ValueError('Unexpected ion-distance trace')
    result = {'species':audit['species'], 'state':'102ps continuation and descriptive diagnostics completed',
              'continuation_duration_ps':100., 'total_duration_ps':102., 'sampled_frames':2001,
              'temperature_mean_K':float(data[:,2].mean()),
              'temperature_second_half_mean_K':float(data[data[:,0]>=52.,2].mean()),
              'bond_range_nm':[min(lengths),max(lengths)], f'{oxygen_label}_oxygen_indices_1based':selected,
              'water_RDF_windows':rdfs,
              f'{oxygen_label}_O_to_Na_minimum_distance_range_nm':[float(distance[:,1].min()),float(distance[:,1].max())],
              'fraction_saved_frames_with_any_O_Na_below_0_35nm':float((distance[:,1]<.35).mean()),
              'contact_scope':'minimum over all selected O/Na sites;0.35nm is a descriptive cutoff,not validated ion-pair equilibrium',
              'concentration_mM':audit['nominal_solute_concentration_mM'],
              'limitations':'one fixed protonation microstate and trajectory per species;correlated windows,not independent replicates;not the experimental concentration or full cPBSzI mixture',
              'physical_approval':False,'full_buffer_ready':False,
              'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [out/'nvt.tpr',out/'nvt.cpt',trajectories[0],out/'source.cpt']}}
    write_json(target/'result.json',result)
    write_json(ROOT/f'models/buffer-parameters/{name.lower()}-solution-102ps-audit.json',result)
    return result


def execute(name, lock_fd):
    source = ROOT/'runs/citrate-refined-controls-v2'/name
    out = BASE/name; out.mkdir(exist_ok=True)
    if (out/'completed.json').exists():
        return
    audit = json.loads((source/'result.json').read_text())
    required = {'md.cpt':'checkpoint_sha256','citrate.itp':'itp_sha256','types.itp':'types_sha256','topol.top':'topology_sha256'}
    if any(sha(source/file)!=audit[key] for file,key in required.items()):
        raise ValueError('Source citrate evidence changed')
    spec = {file:sha(source/file) for file in ['md.tpr','md.cpt','citrate.itp','types.itp','topol.top']}
    if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text()) != spec:
        raise ValueError('Source changed; continuation preserved')
    write_json(out/'spec.json',spec)
    if not (out/'nvt.tpr').exists():
        shutil.copy2(source/'md.cpt',out/'source.cpt')
        run(['convert-tpr','-s',str(source/'md.tpr'),'-extend','100','-o','nvt.tpr'],out,'convert',lock_fd=lock_fd)
    command = ['mdrun','-s','nvt.tpr','-deffnm','nvt','-ntmpi','1','-ntomp','2','-cpt','1','-maxh','1']
    if (out/'nvt.cpt').exists():
        command += ['-cpi','nvt.cpt','-append']
    else:
        if list(out.glob('nvt.part*.log')):
            raise RuntimeError('Partial MD without checkpoint requires review')
        command += ['-cpi','source.cpt','-noappend']
    write_json(BASE/'status.json',{'state':'running','active_species':name,'updated_utc':datetime.now(timezone.utc).isoformat()})
    run(command,out,'launch',lock_fd=lock_fd)
    if any(sha(source/file)!=digest for file,digest in spec.items()):
        raise ValueError('Original files changed during continuation')
    result = analyze(out,source,name,audit)
    write_json(out/'completed.json',result)
    print(name,result['state'],flush=True)


def main():
    BASE.mkdir(exist_ok=True)
    with (BASE/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            for name in ['CIT3','HCIT2']:
                execute(name,lock.fileno())
            write_json(BASE/'status.json',{'state':'completed','physical_approval':False,'full_buffer_ready':False})
        except Exception as error:
            write_json(BASE/'status.json',{'state':'held_preserved','reason':str(error)})
            raise


if __name__ == '__main__':
    main()
