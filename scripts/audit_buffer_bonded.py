"""Typing/bonded-coverage diagnostics only, no charge assignment or MD topology export."""
import argparse,json,os,subprocess,hashlib,re
from rdkit import Chem
from pathlib import Path
from acpype.utils import bundled_amber_dir
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect-existing',action='store_true');args=parser.parse_args()
    out=ROOT/'runs/buffer-bonded-preflight'
    if out.exists() and not args.inspect_existing:raise RuntimeError('Existing diagnostics preserved;use --inspect-existing to review without recomputing')
    if not out.exists() and args.inspect_existing:raise RuntimeError('No existing diagnostics')
    out.mkdir(exist_ok=True);home=Path(bundled_amber_dir());env=os.environ.copy();env['AMBERHOME']=str(home);env['PATH']=str(home/'bin')+':'+env['PATH'];env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    rows=[]
    for item in json.loads((ROOT/'models/buffer-parameters/polyatomic-species-audit.json').read_text())['candidates']:
        name=item['name'];work=out/name;work.mkdir(exist_ok=True);source=ROOT/'models/buffer-parameters'/(name+'.mol');report={'name':name,'formal_charge':item['formal_charge'],'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'scope':'typing and bonded lookup only;no force-field partial charges computed/approved;not an MD topology'}
        commands=[['antechamber','-i',str(source),'-fi','mdl','-o','typed.mol2','-fo','mol2','-at','gaff2','-j','4','-nc',str(item['formal_charge']),'-dr','no','-s','2'],['parmchk2','-i','typed.mol2','-f','mol2','-o','coverage.frcmod','-s','gaff2']]
        try:
            for index,command in enumerate([] if args.inspect_existing else commands):
                with (work/f'command-{index}.log').open('w') as log:subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60,check=True)
            mol2=(work/'typed.mol2').read_text();atoms=mol2.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0];types=[line.split()[5] for line in atoms.splitlines() if line.strip()]
            if len(types)!=item['atoms']:raise ValueError('Typing atom count changed')
            warnings=[line.strip() for log in work.glob('command-*.log') for line in log.read_text().splitlines() if re.search(r'warning|fatal|error:',line,re.I)]
            original=Chem.MolFromMolFile(str(source),removeHs=False)
            table=mol2.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0]
            typed_bonds={tuple(sorted((int(r[1])-1,int(r[2])-1))):r[3] for r in [l.split() for l in table.splitlines() if l.strip()]}
            expected={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(int(b.GetBondTypeAsDouble())) for b in original.GetBonds()}
            if set(typed_bonds)!=set(expected):raise ValueError('Typing changed connectivity')
            changes=[{'atoms_1based':[i+1 for i in edge],'source_order':order,'typed_order':typed_bonds[edge]} for edge,order in expected.items() if typed_bonds[edge]!=order]
            frcmod=(work/'coverage.frcmod').read_text();flags=[line.strip() for line in frcmod.splitlines() if re.search(r'ATTN|need\s+revision',line,re.I)]
            report.update({'state':('blocked:typing warning' if warnings else 'bond-order representation review required' if changes else 'bonded lookup completed'),'typing_warnings':warnings,'bond_order_changes':changes,'atom_types':sorted(set(types)),'missing_or_revision_flags':flags,'analogy_notes':[line.strip() for line in frcmod.splitlines() if 'penalty' in line.lower()],'frcmod_sha256':hashlib.sha256((work/'coverage.frcmod').read_bytes()).hexdigest(),'physical_validation':False})
        except Exception as error:report.update({'state':'diagnostic blocked','reason':str(error)})
        rows.append(report)
    summary={'species':rows,'full_buffer_ready_for_MD':False,'charges_assigned':False,'limitations':'successful lookup does not validate highly charged/polyatomic ions,charge distribution,linear azide geometry,water interactions orpH speciation;diagnostic mol2 charge fields are never exported/used for MD'}
    (ROOT/'models/buffer-parameters/bonded-coverage-audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
