"""Check all saved phosphate bond lengths independently of finite-energy checks."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    path=ROOT/'models/buffer-parameters/projected-phosphate-water-audit.json'
    results=json.loads(path.read_text())
    for record in results:
        out=ROOT/record['run_directory'];n=7 if record['species']=='dihydrogen_phosphate-' else 6
        mol=(out/'phosphate.itp').read_text()
        rows=[s.split(';',1)[0].split() for s in mol.split('[ bonds ]')[1].split('[ pairs ]')[0].splitlines()]
        bonds=[(int(r[0])-1,int(r[1])-1) for r in rows if r]
        target=out/'phosphate-geometry.g96'
        if not target.exists():
            (out/'phosphate.ndx').write_text('[ PHOSPHATE ]\n'+' '.join(map(str,range(1,n+1)))+'\n')
            with (out/'geometry-export.log').open('w') as stream:
                subprocess.run(['gmx','trjconv','-s','md.tpr','-f','md.trr','-n','phosphate.ndx',
                                '-o',target.name,'-ndec','9'],cwd=out,input='PHOSPHATE\n',text=True,
                               stdout=stream,stderr=subprocess.STDOUT,check=True)
        lines=target.read_text().splitlines()
        xyz=np.array([[[float(v) for v in row.split()] for row in lines[i+1:i+n+1]]
                      for i,line in enumerate(lines) if line=='POSITIONRED'])
        assert xyz.shape==(41,n,3) and np.isfinite(xyz).all()
        lengths=[]
        for a,b in bonds:
            delta=xyz[:,a]-xyz[:,b];delta-=3*np.rint(delta/3)
            lengths.extend(np.linalg.norm(delta,axis=1).tolist())
        assert min(lengths)>.07 and max(lengths)<.20
        record['saved_frame_geometry']={'frames':41,'bond_length_range_nm':[min(lengths),max(lengths)],
            'criterion':'all defined OH/PO bonds stay within broad0.07–0.20nm numerical envelope;not an equilibrium bond-distribution validation',
            'trajectory_sha256':hashlib.sha256((out/'md.trr').read_bytes()).hexdigest()}
        (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    path.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps([{'species':r['species'],**r['saved_frame_geometry']} for r in results],indent=2))


if __name__=='__main__':main()
