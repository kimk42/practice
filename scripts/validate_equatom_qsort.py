"""Verify complete path sorting and installed PC4 charge-equivalence behavior."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import resource
import subprocess

import numpy as np
from rdkit import Chem

from pc_length_parameters import ROOT, preflight, sha
from run_pc_sqm_checkpointed import write_json

OUT=ROOT/'runs/equatom-qsort-validation-v1'
BUILD=Path('/workspace/md-tools/equatom-qsort-v1/build.json')


def ac_rows(path):
    rows=[]
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM'):
            continue
        # AC writes adjacent fixed-width atom/residue names: C100MOL is not
        # one atom name. Match the installed writer and pinned upstream reader.
        fields=line[22:].split()
        if len(fields)!=6:
            raise ValueError('Unexpected AC residue/coordinate/charge/type fields')
        number=str(int(line[4:11]));name=line[13:17].strip();residue=line[17:21].strip()
        if not name or not residue:
            raise ValueError('Missing fixed-width AC names')
        rows.append(['ATOM',number,name,residue,*fields])
    return rows


def printed_charges(path,n):
    text=path.read_text().split('Mulliken Charge',1)[1].split('Total Mulliken Charge',1)[0]
    rows=[line.split() for line in text.splitlines() if re.match(r'^\s*\d+\s+[A-Z][a-z]?\s+[-\d.]',line)]
    if len(rows)!=n:
        raise ValueError('Incomplete raw QM charge table')
    return np.array([float(row[-1]) for row in rows])


def average_charges(charges,representatives):
    result=np.empty_like(charges)
    for representative in sorted(set(representatives)):
        indices=[i for i,r in enumerate(representatives) if r==representative]
        # Preserve the upstream atom-index summation order.
        value=sum(float(charges[i]) for i in indices)/len(indices)
        result[indices]=value
    return result


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))


def main():
    OUT.mkdir()
    build=json.loads(BUILD.read_text());results={}
    for dp in [4,8]:
        mol,info=preflight(dp);base=ROOT/f'runs/pc-length-forcefield-v2/pc{dp}'
        atoms=ac_rows(base/'ANTECHAMBER_AC.AC')
        edges=[tuple(sorted(map(int,line.split()[2:4]))) for line in (base/'ANTECHAMBER_AC.AC').read_text().splitlines() if line.startswith('BOND')]
        expected={tuple(sorted([b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1])) for b in mol.GetBonds()}
        if len(atoms)!=mol.GetNumAtoms() or set(edges)!=expected:
            raise ValueError('AC atom count/bond graph differs')
        neighbors=[[] for _ in atoms]
        for i,j in edges:
            neighbors[i-1].append(j-1);neighbors[j-1].append(i-1)
        folder=OUT/f'pc{dp}';folder.mkdir()
        inp=folder/'graph.txt'
        inp.write_text(str(len(atoms))+'\n'+''.join(f'{atom.GetAtomicNum()} {row[2]} '+
            ' '.join(map(str,neighbors[i]+[-1]*(6-len(neighbors[i]))))+'\n'
            for i,(atom,row) in enumerate(zip(mol.GetAtoms(),atoms))))
        outputs={}
        options=[('reference',-1),('qsort',-1)] if dp==4 else [('qsort',-1),('qsort',30)]
        for implementation,length in options:
            binary=Path(build['binaries'][implementation]['path'])
            if sha(binary)!=build['binaries'][implementation]['sha256']:
                raise ValueError('Binary changed')
            label=f'{implementation}-path{length}'
            with (folder/(label+'.log')).open('x') as log:
                subprocess.run([str(binary),str(inp),str(folder/(label+'.json')),str(length)],
                               stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300,preexec_fn=limits)
            outputs[label]=json.loads((folder/(label+'.json')).read_text())
            print(dp,label,'classes',len(set(outputs[label]['representatives_1based'])),
                  'maxpaths',max(outputs[label]['path_counts']),flush=True)
        full=outputs['qsort-path-1']
        charges=printed_charges(base/'accepted-sqm.out',len(atoms))
        averaged=average_charges(charges,full['representatives_1based'])
        ranks=list(Chem.CanonicalRankAtoms(mol,breakTies=False,includeChirality=False))
        expected_representatives=[ranks.index(r)+1 for r in ranks]
        result={'molecule':info,'unlimited_path_classes':len(set(full['representatives_1based'])),
                'max_path_count_per_atom':max(full['path_counts']),
                'canonical_graph_partition_identical':expected_representatives==full['representatives_1based'],
                'raw_QM_sum_e':float(charges.sum()),'averaged_sum_e':float(averaged.sum()),
                'averaged_charges_e':averaged.tolist(),'representatives_1based':full['representatives_1based']}
        if dp==4:
            if outputs['reference-path-1']!=full:
                raise ValueError('Sort replacement changed reference equivalence')
            installed=ac_rows(base/'ANTECHAMBER_AM1BCC_PRE.AC')
            delta=float(np.max(np.abs(averaged-np.array([float(r[8]) for r in installed]))))
            if delta>5.01e-7:
                raise ValueError('Reference differs from installed PC4 averaging')
            result.update(reference_and_qsort_identical=True,installed_PC4_preBCC_max_charge_difference_e=delta)
        else:
            truncated=average_charges(charges,outputs['qsort-path30']['representatives_1based'])
            result['path30_max_charge_difference_e']=float(np.max(np.abs(averaged-truncated)))
            result['path30_partition_identical']=outputs['qsort-path30']['representatives_1based']==full['representatives_1based']
        write_json(folder/'result.json',result);results[str(dp)]=result
    report={'checked_utc':datetime.now(timezone.utc).isoformat(),'build':build,'cases':results,
            'decision':'Use only full-path fast-sort result after PC4 parity validation;do not truncate paths or infer missing QM digits.',
            'physical_forcefield_approval':False}
    write_json(OUT/'result.json',report)
    write_json(ROOT/'models/surface-revision-v2/equatom-qsort-validation-20261030.json',report)


if __name__=='__main__':
    main()
