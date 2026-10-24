"""Plot recorded smoke-test outputs; no simulated or interpolated observations added."""
import csv,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/protein-explorer-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
source=root/'runs/engine-smoke-20261022T093349251988Z/output'
report=json.loads((source/'result.json').read_text())
with (source/'energy.csv').open() as f: rows=list(csv.reader(f))[1:]
time=[int(r[0])*report['config']['timestep_fs']/1000 for r in rows]
energy=[float(r[1]) for r in rows];temp=[float(r[2]) for r in rows]
plt.style.use('seaborn-v0_8-whitegrid')
fig,axes=plt.subplots(1,2,figsize=(12,5))
fig.suptitle('Executed MD: synthetic two-particle engine test',fontsize=17,fontweight='bold',y=.98)
axes[0].plot(time,energy,'o-',color='#137f88');axes[0].set_ylabel('Total potential energy (kJ/mol)');axes[0].set_title('Recorded system energy')
axes[1].plot(time,temp,'o-',color='#aa6832',label='Recorded temperature');axes[1].axhline(310,color='#64748b',ls='--',label='Thermostat target: 310 K');axes[1].set_ylabel('Instantaneous temperature (K)');axes[1].set_title('Short-run temperature fluctuations');axes[1].legend(fontsize=9)
for ax in axes:ax.set_xlabel('Simulation time (ps)');ax.set_xlim(0,1)
fig.text(.5,.055,'CPU | 1,000 steps × 1 fs = 1 ps | 10 recorded frames | seed 42',ha='center',fontsize=10)
fig.text(.5,.012,'No fibrinogen, polymer, or water present. These curves do not measure adsorption or binding strength.',ha='center',fontsize=10,color='#a12d32')
fig.tight_layout(rect=(0,.1,1,.9))
fig.savefig(root/'artifacts/engine-result.png',dpi=160)
fig.savefig(root/'artifacts/engine-result.svg')
import shutil
for name in ('energy.csv','result.json'):shutil.copyfile(source/name,root/'artifacts'/name)
print('Saved artifacts/engine-result.png and .svg, energy.csv, result.json')
