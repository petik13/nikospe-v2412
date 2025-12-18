import numpy as np
import shutil
import subprocess
import helperFuncs as hf
import os

# =============================================================================
logger = hf.setup_logger('log')
hf.console("Starting workflow")
subprocess.run(['./clean.sh'])
shutil.rmtree('0')
shutil.copytree('0.orig', '0')

zeta0 = 0.0207
lam = 7.0
Ucur = 0.0
R0 = 0.5
rampperiod = 3.0
Co = 0.08
Nproc = 32 # Number of processors for parallel run
procD = [8, 4, 1]

steepness = 2 * zeta0/lam
xsponge = 2 * lam
Lsponge = 2 * lam
xbody = 7 * lam 
xdamp = 12 * lam
Lxdamp = 3 * lam 
xmax = 15 * lam 
xmin = 0.0

ydamp = 5.0*lam
Lydamp = 1.5 * lam 
ymin = - 6.5*lam
ymax = 6.5*lam
zmin = -1.0*lam
zmax = 0.0

# Mesh resolution
discX = 20 # cells per lam
discY = discX // 2
discZ = int(1.25*discX)

Nx = int(xmax / lam * discX)
Ny = int((ymax - ymin) / lam * discY)
Nz = int(-zmin / lam * discZ)

# Timestep calculation
k = 2 * np.pi / lam
g = 9.81
omega = np.sqrt(g * k)
period = 2 * np.pi / omega
celerity = omega / k
deltaT = Co * lam / (celerity * discX)
hf.console(f"Calculated timestep deltaT = {deltaT:.6f} s for Co = {Co}")

def update_file(name, value, path='system/blockMeshDict', endl=';'):
    if not os.path.isfile(path):
        logger.warning(f"file not found at {path}, skipping update.")
        return

    with open(path, 'r') as f:
        lines = f.readlines()

    def repl(line: str, key: str, value) -> str:
        stripped = line.lstrip()

        # ignore empty lines and comments
        if not stripped or stripped.startswith('//') or stripped.startswith('/*'):
            return line

        # get first token (keyword)
        first = stripped.split(None, 1)[0]  # e.g. "n" or "numberOfSubdomains"
        if first == key:
            indent = line[:len(line) - len(stripped)]
            return f"{indent}{key}\t{value}{endl}\n"

        return line

    new_lines = [repl(ln, name, value) for ln in lines]

    with open(path, 'w') as f:
        f.writelines(new_lines)

# -- Modify blockMeshDict
hf.console("Modifying blockMeshDict")
update_file('xmax', xmax)
update_file('ymin', ymin)
update_file('ymax', ymax)
update_file('zmin', zmin)
update_file('Nx', Nx)
update_file('Ny', Ny)
update_file('Nz', Nz)
# -- Modify snappyHexMeshDict
hf.console("Modifying snappyHexMeshDict")
update_file('bodyXpos', xbody, path='system/snappyHexMeshDict')

# -- Modify waveCurConditions
wCpath = os.path.join('constant', 'waveCurConditions')
hf.console("Modifying waveCurConditions")
update_file('steepness', steepness, path=wCpath)
update_file('wavelength', lam, path=wCpath)
update_file('Ucur', Ucur, path=wCpath)
update_file('waterdepth', -zmin, path=wCpath)
update_file('xdamp', xdamp, path=wCpath)
update_file('Lxdamp', Lxdamp, path=wCpath)
update_file('ydamp', ydamp, path=wCpath)
update_file('Lydamp', Lydamp, path=wCpath)
update_file('xsponge', xsponge, path=wCpath)
update_file('Lsponge', Lsponge, path=wCpath)
update_file('xc', xbody, path=wCpath)
update_file('R0', R0, path=wCpath)
update_file('rampperiod', rampperiod, path=wCpath)

# -- Modify meanVal.py
mvpath = 'meanVal.py'
hf.console("Modifying meanVal.py")
update_file('lam', f' = {lam}', path=mvpath, endl='')
update_file('zeta', f' = {zeta0}', path=mvpath, endl='')

# -- Modify topoSetDict
tsdpath = os.path.join('system', 'topoSetDict')
hf.console("Modifying topoSetDict")
update_file('xmin', xbody - 2*R0, path=tsdpath)
update_file('xmax', xbody + 2*R0, path=tsdpath)
update_file('ymin', -2*R0, path=tsdpath)
update_file('ymax', 2*R0, path=tsdpath)
update_file('zmin', zmin, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

# -- Modify runCase.sh
rcpath = 'runCase.sh'
hf.console("Modifying runCase.sh")
update_file('Nproc', f'={Nproc}', path=rcpath, endl='')

# -- Modify decomposeParDict
dpdpath = os.path.join('system', 'decomposeParDict')
hf.console("Modifying decomposeParDict")
update_file('numberOfSubdomains', Nproc, path=dpdpath)
update_file('n', f'\t({procD[0]} {procD[1]} {procD[2]}) ', path=dpdpath)

# -- Modify controlDict
cdpath = os.path.join('system', 'controlDict')
hf.console("Modifying controlDict")
update_file('deltaT', f'{deltaT:.5f}', path=cdpath)
update_file('cvPoint', f'({xbody} 0 { -R0})', path=cdpath)
update_file('CofR', f'({xbody} 0 0)', path=cdpath)
# -- Ready to run
subprocess.run(['./mesh.sh'])
subprocess.run(['./runCase.sh'])
