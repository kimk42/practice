"""Compute unapproved AM1-BCC candidates,not measured or validated buffer parameters.
Known azide typing warning is not bypassed;azide is excluded.
"""
import json,subprocess,os,hashlib,concurrent.futures,tempfile
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from acpype.utils import bundled_amber_dir
ROOT=Path(__file__).resolve().parents[1]
def main():
    manifest=ROOT/'runs/buffer-charge-active.json'
    if manifest.exists():raise RuntimeError('Candidate run already registered;inspect/preserve it rather than duplicating')
    out=Path(tempfile.mkdtemp(prefix='buffer-charge-candidates-',dir=ROOT/'runs'));manifest.write_text(json.dumps({'directory':str(out.relative_to(ROOT)),'state':'unapproved candidates running','azide':'excluded due to prior typing warning'},indent=2)+'\n')
    home=Path(bundled_amber_dir());env=os.environ.copy();env['AMBERHOME']=str(home);env['PATH']=str(home/'bin')+':'+env['PATH'];env['LD_LIBRARY_PATH']='/workspace/md-tools/root/usr/lib/x86_64-linux-gnu/openblas-pthread'+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '');env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    species=json.loads((ROOT/'models/buffer-parameters/polyatomic-species-audit.json').read_text())['candidates'];species=[s for s in species if s['name']!='azide-']
    def calculate(item):
        work=out/item['name'];work.mkdir();mol=Chem.AddHs(Chem.MolFromSmiles(item['smiles']));params=AllChem.ETKDGv3();params.randomSeed=2026
        report={'species':item['name'],'expected_formal_charge':item['formal_charge'],'method':'AM1-BCC/GAFF2,isolated gas-phase AM1 optimization,full bond typing with-j4;candidate only;output representation must be reviewed','approved_for_MD':False}
        try:
            if Chem.GetFormalCharge(mol)!=item['formal_charge']:raise ValueError('Input formal charge mismatch')
            if AllChem.EmbedMolecule(mol,params)!=0:raise ValueError('Candidate initial geometry embedding failed')
            Chem.MolToMolFile(mol,str(work/'input.mol'))
            command=['antechamber','-i','input.mol','-fi','mdl','-o','candidate.mol2','-fo','mol2','-c','bcc','-at','gaff2','-j','4','-nc',str(item['formal_charge']),'-m','1','-eq','1','-df','2','-dr','no','-pf','n','-s','2']
            with (work/'calculation.log').open('w') as log:
                process=subprocess.Popen(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:code=process.wait(timeout=1200)
                except subprocess.TimeoutExpired:
                    import signal
                    os.killpg(process.pid,signal.SIGTERM);process.wait();raise
                if code:raise RuntimeError('Antechamber exit '+str(code))
            log=(work/'calculation.log').read_text();warnings=[line for line in log.splitlines() if 'warning' in line.lower()]
            if warnings:raise RuntimeError('Calculation warning requires review:'+str(warnings))
            qm=(work/'sqm.out').read_text()
            if 'Calculation Completed' not in qm:raise RuntimeError('Quantum completion marker absent')
            text=(work/'candidate.mol2').read_text();sections={};section=None
            for line in text.splitlines():
                if line.startswith('@<TRIPOS>'):section=line[9:];sections[section]=[]
                elif section and line.strip():sections[section].append(line.split())
            atoms=sections['ATOM'];bonds=sections['BOND'];numbers={'c':6,'h':1,'o':8,'p':15}
            if [a.GetAtomicNum() for a in mol.GetAtoms()]!=[numbers[r[5][0]] for r in atoms]:raise ValueError('Atom order/element mismatch')
            expected={tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in mol.GetBonds()};actual={tuple(sorted(map(int,r[1:3]))) for r in bonds}
            if expected!=actual:raise ValueError('Connectivity mismatch')
            charge=sum(float(r[8]) for r in atoms)
            if abs(charge-item['formal_charge'])>1e-4:raise ValueError('Net candidate charge mismatch')
            xyz=np.array([[float(x) for x in r[2:5]] for r in atoms]);lengths=[float(np.linalg.norm(xyz[i-1]-xyz[j-1])) for i,j in expected]
            if not np.isfinite(xyz).all() or min(lengths)<.65 or max(lengths)>2.1:raise ValueError('Abnormal candidate bonded geometry')
            report.update({'state':'candidate quantum/charge/connectivity checks passed','net_calculated_charge':charge,'atoms':len(atoms),'bond_length_range_A':[min(lengths),max(lengths)],'molecule_sha256':hashlib.sha256((work/'candidate.mol2').read_bytes()).hexdigest(),'source_input_sha256':hashlib.sha256((work/'input.mol').read_bytes()).hexdigest(),'candidate_directory':str(work.relative_to(ROOT)),'limitations':['gas-phase multiply charged ion electrostatics require independent validation','no hydration/ion-pairing validation or population assignment','no MD topology exported;not part of running baseline or full-buffer system']})
        except Exception as error:report.update({'state':'candidate blocked','reason':str(error),'candidate_directory':str(work.relative_to(ROOT))})
        return report
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(calculate,species))
    report={'state':'unapproved candidate calculations inspected','species':rows,'azide':'not calculated;prior automatic-typing warning remains blocked','full_buffer_ready':False}
    (ROOT/'models/buffer-parameters/charge-candidate-audit.json').write_text(json.dumps(report,indent=2)+'\n');manifest.write_text(json.dumps({'directory':str(out.relative_to(ROOT)),'state':'candidate calculations finished;not approved for MD'},indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
