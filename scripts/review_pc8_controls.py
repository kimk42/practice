"""Read completed PC8 controls independently, without restarting simulations."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from pc_bcc_from_equivalence import validate_equivalence_input
from pc_length_parameters import ROOT, sha
from run_pc_sqm_checkpointed import write_json


def main():
    out=ROOT/'runs/pc8-dry-controls-review-20261030'
    out.mkdir()
    parameter_path=ROOT/'runs/pc-length-forcefield-v2/pc8/parameter-audit.json'
    parameter=json.loads(parameter_path.read_text())
    topology=ROOT/parameter['topology_folder']
    for name,digest in parameter['outputs_sha256'].items():
        if sha(topology/name)!=digest:
            raise ValueError('Reviewed parameter output changed')
    equivalence=json.loads((ROOT/'models/surface-revision-v2/equatom-qsort-validation-20261030.json').read_text())
    source=ROOT/'runs/pc-length-forcefield-v2/pc8'
    validate_equivalence_input(8,equivalence['cases']['8'],source/'accepted-sqm.out',source/'ANTECHAMBER_AC.AC')
    parity=json.loads((ROOT/'runs/pc4-split-bcc-parity-v2/bcc-route-audit.json').read_text())
    if not parity['installed_PC4_final_MOL2_fields_identical']:
        raise ValueError('Installed PC4 BCC parity missing')
    records=[]
    for seed in (2026,2027):
        case=f'pc8-rho0.6-seed{seed}'; folder=ROOT/'runs/pc-bulk-controls-v2'/case
        result=json.loads((folder/'result.json').read_text())
        spec=json.loads((folder/'spec.json').read_text())
        if sha(folder/'pc.itp')!=spec['itp_sha256']:
            raise ValueError('Local topology changed')
        if spec['itp_sha256']!=parameter['outputs_sha256']['PC8_GMX.itp']:
            raise ValueError('Control used a different parameter candidate')
        evidence={}
        for stage,duration,expected_frames in [('nvt',2.,41),('npt',20.,401)]:
            log=(folder/(stage+'.log')).read_text()
            if re.search(r'LINCS WARNING|constraint.*warning|Fatal error',log,re.I):
                raise ValueError('Numerical warning requires review')
            if sha(folder/(stage+'.cpt'))!=result['checkpoints_sha256'][stage+'.cpt']:
                raise ValueError('Completed checkpoint changed')
            xvg=out/(case+'-'+stage+'-energy.xvg')
            selection='Potential\nTemperature\n0\n' if stage=='nvt' else 'Potential\nTemperature\nPressure\nVolume\nDensity\n0\n'
            with (out/(case+'-'+stage+'-energy.log')).open('x') as stream:
                subprocess.run(['gmx','energy','-f',str(folder/(stage+'.edr')),'-o',str(xvg)],
                               input=selection,text=True,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=60)
            data=np.array([[float(x) for x in row.split()] for row in xvg.read_text().splitlines()
                           if row and not row.startswith(('#','@'))])
            if len(data)!=expected_frames or not np.isfinite(data).all() or abs(data[-1,0]-duration)>1e-6:
                raise ValueError('Independent energy/duration check failed')
            item={'duration_ps':float(data[-1,0]),'finite_energy_frames':len(data),
                  'temperature_mean_K':float(data[:,2].mean()),'temperature_final_K':float(data[-1,2])}
            if stage=='npt':
                early=float(data[data[:,0]<=5,5].mean());late=float(data[data[:,0]>=15,5].mean())
                item.update(density_first_5ps_kg_m3=early,density_last_5ps_kg_m3=late,
                            density_window_change_percent=100*(late/early-1))
            evidence[stage]=item
        with (out/(case+'-trajectory-check.log')).open('x') as stream:
            subprocess.run(['gmx','check','-f',str(folder/'npt.xtc')],stdout=stream,
                           stderr=subprocess.STDOUT,check=True,timeout=60)
        text=(out/(case+'-trajectory-check.log')).read_text()
        if not re.search(r'^Coords\s+201\s+0\.1',text,re.M):
            raise ValueError('NPT coordinate frame count differs')
        records.append({'case':case,'EM_Fmax_kJ_mol_nm':result['EM_Fmax_kJ_mol_nm'],
                        'stages':evidence,'NPT_coordinate_frames':201,
                        'result_sha256':sha(folder/'result.json'),
                        'outputs_sha256':{p.name:sha(p) for p in folder.iterdir()
                                          if p.suffix in ('.cpt','.tpr','.edr','.xtc')}})
    preserved=json.loads((ROOT/'runs/pc-length-preservation-before.json').read_text())
    preserved.update(json.loads((ROOT/'runs/pc-continuation-preservation-20261028.json').read_text())['new_checkpoints'])
    for name in ['CIT3','HCIT2']:
        spec=json.loads((ROOT/f'runs/citrate-solution-102ps-v1/{name}/spec.json').read_text())
        preserved[f'runs/citrate-refined-controls-v2/{name}/md.cpt']=spec['md.cpt']
    for name,digest in preserved.items():
        if sha(ROOT/name)!=digest:
            raise ValueError(f'Preserved checkpoint changed: {name}')
    report={'checked_utc':datetime.now(timezone.utc).isoformat(),
            'completed_cases':records,'parameter_audit_sha256':sha(parameter_path),
            'independently_bound_to_original_QM_and_graph':True,
            'installed_PC4_BCC_parity':True,'preserved_prior_checkpoints':len(preserved),
            'preserved_checkpoint_sha256':preserved,
            'scope':'PC8 OH-ended dry numerical controls,initial density0.6g/cm3,300K;not coating/FXII/cPBSzI',
            'interpretation':'Both protocols complete;density drift remains.No equilibrium or physical forcefield approval.',
            'independence':'two packing/velocity seeds;one common intrachain conformer;PC4 initial packing used a different boundary constraint',
            'physical_approval':False,'original_simulations_modified':False}
    write_json(out/'audit.json',report)
    write_json(ROOT/'models/surface-revision-v2/pc8-dry-controls-interim-audit.json',report)
    write_json(ROOT/'models/surface-revision-v2/pc8-recovered-parameter-audit.json',parameter)
    print(json.dumps({'completed_cases':len(records),'preserved_checkpoints':len(preserved),
                      'density_window_change_percent':[r['stages']['npt']['density_window_change_percent'] for r in records]}))


if __name__=='__main__':
    main()
