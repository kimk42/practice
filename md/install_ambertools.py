"""Install checksum-pinned AmberTools locally without registering in the user home."""
from pathlib import Path
import argparse
import hashlib,io,os,subprocess,sys,tempfile,urllib.request,zipfile,tarfile
import zstandard
ROOT=Path(__file__).resolve().parents[1]
TOOLS=Path(os.environ.get('FIBRINOGEN_TOOLS_DIR',ROOT/'.tools')).expanduser().resolve()
PACKAGE='https://conda.anaconda.org/conda-forge/linux-64/conda-standalone-26.7.0-h643be8f_0.conda'
EXPECTED='1bdee371f8d53752b1af49aa1dcf78f11cb7217c4e6903e2e40dcf7db766f91a'
try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def install(prefix=None):
    prefix=Path(prefix or amber_home()).expanduser().resolve()
    if (prefix/'bin/antechamber').exists() and list((prefix/'conda-meta').glob('ambertools-24.8-*.json')):
        print('AmberTools 24.8 already installed:',prefix);return
    installer=TOOLS/'conda-standalone/standalone_conda/conda.exe'
    # Revalidate the archive even when executable extraction is reusable.
    archive=TOOLS/'installers/conda-standalone-26.7.0-h643be8f_0.conda';archive.parent.mkdir(parents=True,exist_ok=True)
    if not archive.exists():
        with urllib.request.urlopen(PACKAGE,timeout=120) as response:archive.write_bytes(response.read())
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=EXPECTED:raise ValueError('Conda standalone checksum failed')
    with zipfile.ZipFile(archive) as z:
        name=next(n for n in z.namelist() if n.startswith('pkg-') and n.endswith('.tar.zst'))
        with z.open(name) as compressed,zstandard.ZstdDecompressor().stream_reader(compressed) as reader,tarfile.open(fileobj=reader,mode='r|') as tar:
            for member in tar:
                if member.name=='standalone_conda/conda.exe':
                    installer.parent.mkdir(parents=True,exist_ok=True)
                    with installer.open('wb') as f:
                        import shutil
                        shutil.copyfileobj(tar.extractfile(member),f)
                    installer.chmod(0o755);break
            else:raise ValueError('Missing executable in trusted conda package')
    env={**os.environ,'CONDA_REGISTER_ENVS':'false','CONDA_SOLVER':'classic','CONDA_PKGS_DIRS':str(TOOLS/'pkgs'),'CONDA_ENVS_PATH':str(TOOLS/'envs'),'XDG_CACHE_HOME':str(TOOLS/'cache')}
    subprocess.run([str(installer),'create','--yes','--prefix',str(prefix),'--file',str(ROOT/'environment/ambertools-linux64.explicit')],env=env,check=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prefix');args=parser.parse_args();install(args.prefix)
