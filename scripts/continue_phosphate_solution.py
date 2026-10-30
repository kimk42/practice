"""Continue unchanged projected-phosphate controls for descriptive water checks."""
from datetime import datetime, timezone
import fcntl
import json
import shutil

from rdkit import Chem

from pc_length_parameters import ROOT, sha
from run_pc_sqm_checkpointed import write_json
from continue_citrate_solution import run, analyze

BASE=ROOT/'runs/phosphate-solution-102ps-v1'


def execute(name,lock_fd):
    source=ROOT/f'runs/{name.lower()}-projected-water-control'
    out=BASE/name;out.mkdir(exist_ok=True)
    if (out/'completed.json').exists():
        return
    audit=json.loads((source/'result.json').read_text())
    for file,key in [('md.cpt','checkpoint_sha256'),('phosphate.itp','itp_sha256'),('types.itp','types_sha256')]:
        if sha(source/file)!=audit[key]:
            raise ValueError('Phosphate source evidence changed')
    spec={file:sha(source/file) for file in ['md.tpr','md.cpt','phosphate.itp','types.itp','topol.top']}
    if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:
        raise ValueError('Phosphate source changed;continuation preserved')
    write_json(out/'spec.json',spec)
    if not (out/'nvt.tpr').exists():
        shutil.copy2(source/'md.cpt',out/'source.cpt')
        run(['convert-tpr','-s',str(source/'md.tpr'),'-extend','100','-o','nvt.tpr'],out,'convert',lock_fd=lock_fd)
    command=['mdrun','-s','nvt.tpr','-deffnm','nvt','-ntmpi','1','-ntomp','2','-cpt','1','-maxh','1']
    if (out/'nvt.cpt').exists():
        command+=['-cpi','nvt.cpt','-append']
    else:
        if list(out.glob('nvt.part*.log')):
            raise RuntimeError('Partial MD without checkpoint requires review')
        command+=['-cpi','source.cpt','-noappend']
    write_json(BASE/'status.json',{'state':'running','active_species':name,'updated_utc':datetime.now(timezone.utc).isoformat()})
    run(command,out,'launch',lock_fd=lock_fd)
    if any(sha(source/file)!=digest for file,digest in spec.items()):
        raise ValueError('Source changed during continuation')
    mol=Chem.MolFromMolFile(str(ROOT/f'models/buffer-parameters/{audit["species"]}.mol'),removeHs=False)
    # Include every nonprotonated P-bound oxygen, irrespective of the equivalent
    # resonance drawing; exclude hydroxyl oxygen carrying an explicit H.
    oxygens=[a.GetIdx()+1 for a in mol.GetAtoms() if a.GetAtomicNum()==8
             and any(n.GetAtomicNum()==15 for n in a.GetNeighbors())
             and not any(n.GetAtomicNum()==1 for n in a.GetNeighbors())]
    if len(oxygens)!=(2 if name=='PH1' else 3):
        raise ValueError('Unexpected phosphate protonation graph')
    audit['nominal_solute_concentration_mM']=1/6.02214076e23/(27e-24)*1000
    result=analyze(out,source,name,audit,selected_oxygens=oxygens,oxygen_label='nonprotonated_phosphate',maximum_bond_nm=.20)
    result['charge_model_scope']='unchanged previously projected phosphate charges;not a new charge fit'
    write_json(out/'completed.json',result)
    write_json(ROOT/f'models/buffer-parameters/{name.lower()}-solution-102ps-audit.json',result)
    print(name,result['state'],flush=True)


def main():
    BASE.mkdir(exist_ok=True)
    with (BASE/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            for name in ['PH1','PH2']:
                execute(name,lock.fileno())
            write_json(BASE/'status.json',{'state':'completed','physical_approval':False,'full_buffer_ready':False})
        except Exception as error:
            write_json(BASE/'status.json',{'state':'held_preserved','reason':str(error)})
            raise


if __name__=='__main__':
    main()
