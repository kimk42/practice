"""Extract the last complete periodic EM coordinate frame, if one exists."""
import re,subprocess
from pathlib import Path
def last_em_frame(directory:Path,prefix='em'):
    trajectory=directory/(prefix+'.trr');tpr=directory/(prefix+'.tpr')
    if not trajectory.exists() or not tpr.exists():return None
    result=subprocess.run(['gmx','check','-f',str(trajectory)],capture_output=True,text=True)
    matches=re.findall(r'Last frame\s+\d+\s+time\s+([\d.eE+-]+)',result.stderr+result.stdout)
    if not matches:return None
    output=directory/'recovered-em-coordinate.gro'
    with (directory/'recovery-extract.log').open('w') as log:
        subprocess.run(['gmx','trjconv','-f',str(trajectory),'-s',str(tpr),'-dump',matches[-1],'-o',str(output)],input='0\n',text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    return output
