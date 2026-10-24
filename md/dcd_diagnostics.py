"""Read the full-coordinate, orthorhombic OpenMM DCD used by this project.
Coordinates in DCD are Å. Internal diagnostic coordinates are nm.
"""
import csv,struct
from pathlib import Path
import numpy as np
try:
    from .protein_pilot import aligned_rmsd
except ImportError:
    from protein_pilot import aligned_rmsd


def read_record(f):
    raw=f.read(4)
    if len(raw)!=4:raise ValueError('Missing DCD record')
    n=struct.unpack('<i',raw)[0]
    if not 0<n<100000000:raise ValueError('Unexpected DCD record size')
    data=f.read(n);end=f.read(4)
    if len(data)!=n or len(end)!=4 or struct.unpack('<i',end)[0]!=n:raise ValueError('Truncated/malformed DCD record')
    return data


def frames(path):
    with Path(path).open('rb') as f:
        header=read_record(f)
        if header[:4]!=b'CORD' or len(header)!=84:raise ValueError('Not a supported OpenMM coordinate DCD')
        count,start,interval=struct.unpack('<3i',header[4:16])
        has_cell=struct.unpack('<i',header[44:48])[0]
        if not has_cell:raise ValueError('Periodic unit cell required')
        read_record(f);n=struct.unpack('<i',read_record(f))[0]
        for index in range(count):
            cell=np.frombuffer(read_record(f),dtype='<f8')
            if len(cell)!=6 or not np.allclose(cell[[1,3,4]],0,atol=1e-6):raise ValueError('Only orthorhombic DCD cells supported')
            box=cell[[0,2,5]]/10
            axes=[np.frombuffer(read_record(f),dtype='<f4') for _ in range(3)]
            if any(len(a)!=n for a in axes):raise ValueError('DCD particle count mismatch')
            yield start+index*interval,np.column_stack(axes)/10,box


def trajectory_checks(path,prm,reference,output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    atoms=list(prm.topology.atoms());ca=np.array([a.index for a in atoms if a.residue.index<252 and a.name=='CA'])
    chains=np.array([[a.index for a in r.atoms() if a.element.symbol=='C'] for r in prm.topology.residues() if r.name=='PEO'])
    values=[]
    for step,p,box in frames(path):
        if len(p)!=43288 or not np.isfinite(p).all():raise ValueError('Unexpected/non-finite DCD coordinates')
        c=p[ca];rg=float(np.sqrt(np.mean(np.sum((c-c.mean(axis=0))**2,axis=1))))
        d=p[chains[:,-1]]-p[chains[:,0]];d-=box*np.round(d/box);length=np.linalg.norm(d,axis=1);u=d/length[:,None]
        order=float(np.linalg.eigvalsh(1.5*np.einsum('ni,nj->ij',u,u)/len(u)-.5*np.eye(3)).max())
        values.append([step,step*.002,rg,aligned_rmsd(c,reference[ca])*10,float(length.mean()),float(length.min()),order])
    with (out/'trajectory-structure.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['step','extension_time_ps','CA_radius_of_gyration_nm','CA_RMSD_from_extension_start_A','PE_mean_end_to_end_nm','PE_min_end_to_end_nm','PE_nematic_order_max_eigenvalue']);w.writerows(values)
    return {'frames_read':len(values),'last_step':values[-1][0],'CA_Rg_range_nm':[min(x[2] for x in values),max(x[2] for x in values)],'PE_alignment_range':[min(x[-1] for x in values),max(x[-1] for x in values)],'note':'Coarse structural diagnostics from actual DCD frames; not native-state, crystallinity or equilibrium validation.'}
