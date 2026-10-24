"""Inspect uploaded 3GHG metadata and export observed assembly-1 protein C-alpha coordinates for viewing only."""
import json,hashlib,os
from pathlib import Path
from collections import Counter
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
ROOT=Path(__file__).resolve().parents[1]
def main():
    source=ROOT/'data/3ghg-upload/source.cif'
    c=MMCIF2Dict(str(source))
    if c['_entry.id'][0]!='3GHG':raise ValueError('Expected 3GHG')
    idx=c['_pdbx_struct_assembly_gen.assembly_id'].index('1')
    if c['_pdbx_struct_assembly_gen.oper_expression'][idx]!='1':raise ValueError('Nonidentity operation requires explicit assembly reconstruction')
    for i in range(1,4):
        for j in range(1,4):
            if float(c[f'_pdbx_struct_oper_list.matrix[{i}][{j}]'][0])!=float(i==j):raise ValueError('Nonidentity matrix')
        if float(c[f'_pdbx_struct_oper_list.vector[{i}]'][0])!=0:raise ValueError('Nonzero transform')
    selected=set(c['_pdbx_struct_assembly_gen.asym_id_list'][idx].split(','))
    entity=dict(zip(c['_struct_asym.id'],c['_struct_asym.entity_id']))
    chains=[x for x in c['_struct_asym.id'] if x in selected and entity[x] in ('1','2','3')]
    expected=Counter(c['_entity_poly_seq.entity_id'])
    rows=[];seen=set()
    for i,name in enumerate(c['_atom_site.label_atom_id']):
        chain=c['_atom_site.label_asym_id'][i];seq=c['_atom_site.label_seq_id'][i]
        if chain not in chains or name!='CA' or c['_atom_site.pdbx_PDB_model_num'][i]!='1':continue
        if (chain,seq) in seen:continue
        seen.add((chain,seq))
        rows.append({'chain':chain,'seq':int(seq),'auth_seq':int(c['_atom_site.auth_seq_id'][i]),'name':c['_atom_site.label_comp_id'][i],
                     'xyz':[float(c[f'_atom_site.Cartn_{k}'][i]) for k in 'xyz']})
    summary={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'title':c['_struct.title'][0],
             'resolution_A':float(c['_refine.ls_d_res_high'][0]),'assembly_id':'1','representation':'observed C-alpha only; not MD input',
             'chains':[{'chain':ch,'entity':entity[ch],'sequence_residues':expected[entity[ch]],'observed_CA':sum(r['chain']==ch for r in rows),
                        'unobserved_sequence_positions':sorted(set(range(1,expected[entity[ch]]+1))-{r['seq'] for r in rows if r['chain']==ch})} for ch in chains]}
    output=ROOT/'data/3ghg-upload';(output/'assembly1-summary.json').write_text(json.dumps(summary,indent=2))
    lines=['REMARK 900 OBSERVED CA ONLY. VIEWING EXPORT; NOT A SIMULATION MODEL.']
    for i,r in enumerate(rows,1):
        x,y,z=r['xyz'];lines.append(f"ATOM  {i:5d}  CA  {r['name']:3s} {r['chain']}{r['auth_seq']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    (output/'assembly1-view.pdb').write_text('\n'.join(lines)+'\nEND\n')
    os.environ.setdefault('MPLCONFIGDIR','/tmp/fibrinogen-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig=plt.figure(figsize=(12,6));ax=fig.add_subplot(111,projection='3d')
    colors=['#d88534','#4a9cb3','#9362b5','#dfb85b','#36b4a0','#da789d']
    for ch,color in zip(chains,colors):
        part=[r for r in rows if r['chain']==ch];xyz=np.array([r['xyz'] for r in part])
        ax.scatter(*xyz.T,s=2,color=color)
        for a,b in zip(part,part[1:]):
            if b['seq']==a['seq']+1 and np.linalg.norm(np.array(a['xyz'])-b['xyz'])<5:
                ax.plot(*np.array([a['xyz'],b['xyz']]).T,lw=.55,color=color)
        ax.plot([],[],[],color=color,label=f"{ch} ({ {'1':'alpha','2':'beta','3':'gamma'}[entity[ch]]})")
    xyz=np.array([r['xyz'] for r in rows]);ax.set_box_aspect(np.maximum(np.ptp(xyz,axis=0),1))
    ax.set_xlabel('X (Å)');ax.set_ylabel('Y (Å)');ax.set_zlabel('Z (Å)');ax.legend(loc='upper left',fontsize=8)
    ax.view_init(elev=23,azim=-65)
    fig.suptitle('Uploaded 3GHG: observed fibrinogen structure',fontsize=17)
    fig.text(.5,.065,'Assembly 1 | six protein chains | observed C-alpha coordinates only',ha='center')
    fig.text(.5,.025,'Missing regions are not reconstructed. Knob peptides, glycans and calcium are omitted from this VIEW only. No adsorption simulation yet.',ha='center',fontsize=9,color='#a12d32')
    fig.savefig(ROOT/'artifacts/fibrinogen-3ghg.png',dpi=170)
    print(json.dumps({**summary,'chains':[{k:v for k,v in r.items() if k!='unobserved_sequence_positions'} for r in summary['chains']]},indent=2))
if __name__=='__main__':main()
