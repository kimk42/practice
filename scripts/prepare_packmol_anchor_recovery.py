"""Prepare new attempts for unfinished PC cases using the validated PBC variant."""
import fcntl
import json
from pathlib import Path
import re
import shutil
from pc_length_parameters import ROOT,sha
from run_pc_sqm_checkpointed import write_json
from run_pc_packing_margin import BASE,PLAN


def main():
    verified=json.loads((ROOT/'models/surface-revision-v2/packmol-anchor-validation.json').read_text())
    build=json.loads((ROOT/'models/surface-revision-v2/packmol-anchor-build.json').read_text())
    if verified['binary_sha256']!=build['binary_sha256'] or sha(Path(build['binary']))!=build['binary_sha256']:
        raise ValueError('Variant build/validation differs')
    records=[]
    with (BASE/'queue.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for case in json.loads(PLAN.read_text())['cases']:
            root=BASE/case['name'];state=json.loads((root/'status.json').read_text())
            if state['state']=='geometry_passed':continue
            if state['state'] not in ('blocked_preserved','interrupted_preserved','prepared'):
                raise RuntimeError('Unreviewed or active case: '+case['name'])
            with (root/'case.lock').open('a') as case_lock:
                fcntl.flock(case_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                source=root/state.get('attempt_directory','.')
                out=root/'anchor-0001';out.mkdir()
                write_json(out/'previous-status.json',state)
                shutil.copy2(source/'molecule.pdb',out/'molecule.pdb')
                text=re.sub(r'^restart_from .*\n','',(source/'packmol.inp').read_text(),flags=re.M)
                restart=source/'packing.restart'
                if not restart.exists():restart=source/'source.restart'
                provenance=None
                if restart.exists():
                    if len(restart.read_text().splitlines())!=case['chains']:raise ValueError('Incomplete saved poses')
                    shutil.copy2(restart,out/'source.restart');text='restart_from source.restart\n'+text
                    provenance={'source':str(restart.relative_to(ROOT)),'sha256':sha(restart),'scope':'continue last saved poses;no coordinate reset'}
                elif state['state']!='prepared':raise RuntimeError('Attempted case has no saved poses')
                (out/'packmol.inp').write_text(text)
                updated=dict(case,input_sha256=sha(out/'packmol.inp'),restart=provenance,
                    generator_binary=build['binary'],generator_sha256=build['binary_sha256'],
                    boundary='PBC interatomic distances retained;automatic box restriction only on atom1 per molecule')
                write_json(out/'case.json',updated)
                priority=0 if case['DP']==16 and case['seed']==2027 and case['starting_density_g_cm3']==.6 else (
                    1 if case['DP']==8 and case['seed']==2027 and case['starting_density_g_cm3']==.6 else 10)
                write_json(root/'status.json',{'state':'prepared_anchor_variant','attempt_directory':out.name,
                    'recovery_method':'validated_anchor_only_PBC_generator','recovery_priority':priority,'MD_started':False})
                records.append({'case':case['name'],'attempt_directory':str(out.relative_to(ROOT)),
                                'previous_state':state['state'],'restart':provenance})
    write_json(ROOT/'models/surface-revision-v2/packmol-anchor-recovery-plan.json',{
        'cases':records,'binary_sha256':build['binary_sha256'],'completed_PC4_cases_preserved':2,
        'independent_distance_gate_A':2.,'PDB_rounding_allowance_A':.002,'physical_approval':False,
        'supersedes':'v2-v4 claim that selecting atom1 in the input removed whole-molecule containment;upstream automatic restriction persisted'})
    print('Prepared',len(records),'new attempts;completed cases preserved')


if __name__=='__main__':main()
