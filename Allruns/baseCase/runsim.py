import numpy as np
import shutil
import subprocess
import helperFuncs as hf
import os
import argparse
import helpers

# flags:
meshFlag = True
runFlag = True
postFlag = True
# =============================================================================
logger = hf.setup_logger('log')
hf.console("Starting workflow")
subprocess.run(['./clean.sh'])
shutil.rmtree('0')
shutil.copytree('0.orig', '0')
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--lam", type=float, required=True)
parser.add_argument("--heading", type=float, required=True)
parser.add_argument("--Ucur", type=float, required=True)
args = parser.parse_args()

# -- Define parameters
steepness = 1/60
heading = args.heading
lam = args.lam
zeta0 = 0.5*steepness*lam
L_2 = 2.721
B = 51/63.65
depth = 14.5/63.65 # = 0.228
Ucur = args.Ucur
R0 = 0.5
rampperiod = 3.0
Co = 0.2
Nproc = 48 # Number of processors for parallel run
procD = [8, 6, 1]

head_ang = heading*np.pi/180
xsponge = 2 * lam
Lsponge = 2 * lam
xbody = 6 * lam 
xdamp = 12 * lam
Lxdamp = 3 * lam 
xmax = 15 * lam 
xmin = 0.0

ydamp = 7*lam
Lydamp = 2 * lam 
ymin = - 9*lam
ymax = 9*lam
zmin = min(-0.5*lam, -7*depth)
zmax = 0.0

# Mesh resolution
discX = 10 # cells per lam
discY = discX // 2.0
discZ = int(2.0*discX)
Nref = 3 # number of refinements for snappyHexMesh in x-direction.
if lam >= 6.0:
    Nref = 4

Nx = int(xmax / lam * discX)
Ny = int((ymax - ymin) / lam * discY)
Nz = int(-zmin / lam * discZ)

# Timestep calculation
k = 2 * np.pi / lam
g = 9.81
omega = np.sqrt(g * k)
period = 2 * np.pi / omega
celerity = omega / k
Cgroup = 0.5*celerity
T = 2*np.pi/omega
endTime = 3*T + 1.1*(xbody + xbody)/Cgroup
deltaT = Co * lam / ((celerity+Ucur) * discX)/(Nref + 1)
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
update_file('head_ang', -head_ang, path=wCpath)
update_file('steepness', steepness, path=wCpath)
update_file('wavelength', lam, path=wCpath)
update_file('currentspeed', Ucur, path=wCpath)
update_file('waterdepth', -zmin, path=wCpath)
update_file('xdamp', xdamp, path=wCpath)
update_file('Lxdamp', Lxdamp, path=wCpath)
update_file('ydamp', ydamp, path=wCpath)
update_file('Lydamp', Lydamp, path=wCpath)
update_file('xsponge', xsponge, path=wCpath)
update_file('Lsponge', Lsponge, path=wCpath)
update_file('xc', xbody/lam, path=wCpath)
update_file('R0', R0, path=wCpath)
update_file('rampperiod', rampperiod, path=wCpath)

# -- Modify meanVal.py
mvpath = 'meanVal.py'
hf.console("Modifying meanVal.py")
update_file('lam', f' = {lam}', path=mvpath, endl='')
update_file('zeta', f' = {zeta0}', path=mvpath, endl='')
update_file('heading', f' = {heading}', path=mvpath, endl='')
update_file('Ucur', f' = {Ucur}', path=mvpath, endl='')


# -- Modify topoSetDict
tsdpath = os.path.join('system', 'topoSetDict')
hf.console("Modifying topoSetDict")
l1 = 2.0
update_file('xmin', xbody - l1*L_2, path=tsdpath)
update_file('xmax', xbody + l1*L_2, path=tsdpath)
update_file('ymin', -l1*L_2, path=tsdpath)
update_file('ymax', l1*L_2, path=tsdpath)
update_file('zmin', zmin, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

# -- Modify topoSetDict_2
tsdpath = os.path.join('system', 'topoSetDict_2')
hf.console("Modifying topoSetDict_2")
l1 = 3.0
update_file('xmin', xbody - l1*L_2, path=tsdpath)
update_file('xmax', xbody + l1*L_2, path=tsdpath)
update_file('ymin', -l1*L_2, path=tsdpath)
update_file('ymax', l1*L_2, path=tsdpath)
update_file('zmin', zmin, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

# -- Modify topoSetDict.1
update_file('ybox', 0.8*ymax, path='system/topoSetDict.1')
update_file('boxstart', 0, path='system/topoSetDict.1')
update_file('boxend', xdamp, path='system/topoSetDict.1')
update_file('zbox', max(2.5*depth, 0.25*lam), path='system/topoSetDict.1')
# -- Modify topoSetDict.2
update_file('ybox', 0.7*ymax, path='system/topoSetDict.2')
update_file('boxstart', 0, path='system/topoSetDict.2')
update_file('boxend', xdamp, path='system/topoSetDict.2')
update_file('zbox', max(2.0*depth, 0.2*lam), path='system/topoSetDict.2')

# -- Modify topoSetDict.3
l4 = 3.0
update_file('ybox', l4*L_2, path='system/topoSetDict.3')
update_file('boxstart', xbody - l4*L_2, path='system/topoSetDict.3')
update_file('boxend', xbody + l4*L_2, path='system/topoSetDict.3')
update_file('zbox', 1.6*depth, path='system/topoSetDict.3')


# -- Modify topoSetDict.4
l2 = 2.6
p1 = np.array([- l2/2*L_2, l2/2*B, 0])
p2 = np.array([l2/2*L_2, l2/2*B, 0])
p3 = np.array([- l2/2*L_2, -l2/2*B, 0])

rot_ang = -head_ang
rot_mat = np.array([[np.cos(rot_ang), -np.sin(rot_ang)], [np.sin(rot_ang), np.cos(rot_ang)]])
p1_rot = rot_mat @ p1[:2]
p2_rot = rot_mat @ p2[:2]
p3_rot = rot_mat @ p3[:2]

p1_rot[0] += xbody
p2_rot[0] += xbody
p3_rot[0] += xbody

i_ = p2_rot - p1_rot
j_ = p3_rot - p1_rot
k_z = - 1.4*depth

update_file('xbox', p1_rot[0], path='system/topoSetDict.4')
update_file('ybox', p1_rot[1], path='system/topoSetDict.4')
update_file('zbox', 0.0, path='system/topoSetDict.4')
update_file('i_x', i_[0], path='system/topoSetDict.4')
update_file('i_y', i_[1], path='system/topoSetDict.4')
update_file('j_x', j_[0], path='system/topoSetDict.4')
update_file('j_y', j_[1], path='system/topoSetDict.4')
update_file('k_z', k_z, path='system/topoSetDict.4')


# -- Modify topoSetDict.5
l2 = 2.3
p1 = np.array([- l2/2*L_2, l2/2*B, 0])
p2 = np.array([l2/2*L_2, l2/2*B, 0])
p3 = np.array([- l2/2*L_2, -l2/2*B, 0])

rot_ang = -head_ang
rot_mat = np.array([[np.cos(rot_ang), -np.sin(rot_ang)], [np.sin(rot_ang), np.cos(rot_ang)]])
p1_rot = rot_mat @ p1[:2]
p2_rot = rot_mat @ p2[:2]
p3_rot = rot_mat @ p3[:2]

p1_rot[0] += xbody
p2_rot[0] += xbody
p3_rot[0] += xbody

i_ = p2_rot - p1_rot
j_ = p3_rot - p1_rot
k_z = - 1.3*depth

update_file('xbox', p1_rot[0], path='system/topoSetDict.5')
update_file('ybox', p1_rot[1], path='system/topoSetDict.5')
update_file('zbox', 0.0, path='system/topoSetDict.5')
update_file('i_x', i_[0], path='system/topoSetDict.5')
update_file('i_y', i_[1], path='system/topoSetDict.5')
update_file('j_x', j_[0], path='system/topoSetDict.5')
update_file('j_y', j_[1], path='system/topoSetDict.5')
update_file('k_z', k_z, path='system/topoSetDict.5')


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
update_file('endTime', f'{endTime:.2f}', path=cdpath)
update_file('writeInterval', f'{endTime - 1:.2f}', path=cdpath)
update_file('cvPoint', f'({xbody} 0 { -R0})', path=cdpath)
update_file('CofR', f'({xbody} 0 0)', path=cdpath)

# -- Modify rigidBodyMotionProperties
rgpath = os.path.join('constant', 'rigidBodyMotionProperties')
update_file('xG', f'{xbody:.4}', path = rgpath)
update_file('heading', f'{np.pi-head_ang:.4}', path = rgpath)

# -- Modify linMotions.py
rgpath = 'linMotions.py'
update_file('lam', f' = {lam:.2}', path = rgpath,  endl='')
update_file('zeta', f' = {zeta0:.4}', path = rgpath,  endl='')
update_file('heading', f' = {heading:.4}', path = rgpath,  endl='')
update_file('Ucur', f' = {Ucur:.4}', path = rgpath,  endl='')

# -- Modify 0.orig/PhiCur
phicurpath = os.path.join('0.orig', 'PhiCur')
update_file('Ucur', f'{Ucur:.4}', path = phicurpath)
update_file('heading', f'{-heading:.4}', path = phicurpath)

# -- translate surface
subprocess.run(['surfaceTransformPoints', '-rotate', f'((-1 0 0) ({-np.cos(head_ang)} {np.sin(head_ang)} 0))', 'constant/triSurface/DTC.stl', 'constant/triSurface/DTC_rotated.stl'])
subprocess.run(['surfaceTransformPoints', '-translate', f'({xbody} 0 0)', 'constant/triSurface/DTC_rotated.stl', 'constant/triSurface/DTC_moved.stl'])

# -- Ready to run
if meshFlag:
    helpers.mesh(lam)
if runFlag:
    rc = helpers.run_case(Nproc)
    print("shipFlow finished rc=", rc)
if postFlag:
    helpers.post_process()

