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
# Control-volume placement (sphere_diag tests A-G, 2026-10-15): the midfield drift
# force is converged at the base resolution only on a control surface ~8 cells or
# more from the hull and entirely inside the finest refinement box, also in depth.
# 0.1 L from the hull (3 cells) over-predicted the surge drift by 8 %.
parser.add_argument("--cvClear", type=float, default=0.25,
                    help="clearance of the control volume around the hull, sides and below the keel [m]")
parser.add_argument("--zMargin", type=float, default=0.1,
                    help="finest refinement boxes at least this much deeper than the control volume [m]")
args = parser.parse_args()

# -- Define parameters
hull = 'sphere'
steepness = 1/50
heading = args.heading
lam = args.lam
zeta0 = 0.5*steepness*lam
scale = 100.0
L_2 = 0.5
B = 1.0
draft = 0.5
Ucur = args.Ucur
rampperiod = 3.0
Co = 0.2
Nproc = 56 # Number of processors for parallel run
procD = [8, 6, 1]

beta = 0.0

head_ang = heading*np.pi/180
xsponge = 2 * lam
Lsponge = 2 * lam
xbody = 7 * lam 
xdamp = 12 * lam
Lxdamp = 2 * lam 
xmax = 14 * lam 
xmin = 0.0

ydamp = 5*lam
Lydamp = 2 * lam 
ymin = - 7*lam
ymax = 7*lam
zmin = min(-0.5*lam, -7*draft)
zmax = 0.0

# Mesh resolution
discX = 7 # cells per lam
discY = discX // 1.0
discZ = int(2.0*discX)
Nref = 2 # number of refinements for snappyHexMesh in x-direction.
if lam >= 1.0:
    Nref = 3
if lam >= 2.0:
    Nref = 4
if lam >= 2.75:
    Nref = 5
if lam >= 3.5:
    Nref = 6

k = 2 * np.pi / lam
g = 9.81
we  = np.sqrt(g*k*np.tanh(k*(-zmin))) + k*Ucur*np.cos(head_ang) # encounter frequency
tau = we * Ucur / g                       # signed; we = omega_encounter
if tau > 0:
    ratio = ((np.sqrt(1+4*tau) - 1) / (1 - np.sqrt(1-4*tau)))**2
else:
    ratio = 1.0                           # opposing current: scattered waves are longer
lam_mesh = lam

Nx = int(xmax / lam_mesh * discX)
Ny = int((ymax - ymin) / lam_mesh * discY)
Nz = int(-zmin / lam_mesh * discZ)

# Timestep calculation
k = 2 * np.pi / lam
g = 9.81
omega = np.sqrt(g * k)
period = 2 * np.pi / omega
celerity = omega / k
Cgroup = 0.5*celerity
T = 2*np.pi/omega
endTime = rampperiod*T + 0.5*(xbody + xbody)/Cgroup
deltaT = Co * lam_mesh / ((celerity+Ucur) * discX)/(Nref + 1)
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
wCpath = os.path.join('constant', 'waveConditions')
hf.console("Modifying waveConditions")
update_file('headingAngle', -head_ang, path=wCpath)
update_file('steepness', steepness, path=wCpath)
update_file('waveLength', lam, path=wCpath)
update_file('currentSpeed', Ucur, path=wCpath)
update_file('waterDepth', -zmin, path=wCpath)
update_file('xOutlet', xdamp, path=wCpath)
update_file('LOutlet', Lxdamp, path=wCpath)
update_file('ySide', ydamp, path=wCpath)
update_file('LSide', Lydamp, path=wCpath)
update_file('xInlet', xsponge, path=wCpath)
update_file('LInlet', Lsponge, path=wCpath)
update_file('rampPeriods', rampperiod, path=wCpath)
update_file('swaySpeed', np.tan(np.radians(beta)) * Ucur, path=wCpath)

# -- Modify meanVal.py
mvpath = 'meanVal.py'
hf.console("Modifying meanVal.py")
update_file('lam', f' = {lam}', path=mvpath, endl='')
update_file('zeta', f' = {zeta0}', path=mvpath, endl='')
update_file('heading', f' = {heading}', path=mvpath, endl='')
update_file('Ucur', f' = {Ucur}', path=mvpath, endl='')

####### CONTROL VOLUMES ########
# -- Modify topoSetDict
tsdpath = os.path.join('system', 'topoSetDict')
hf.console("Modifying topoSetDict")
# Box with its sides along the mesh axes, around the rotated hull: the L x B
# waterplane rectangle about midship (bow along (-cos h, sin h)), plus a
# clearance of cv_clear*L on each side
# clearance args.cvClear [m] on each side and below the keel (was 0.1 L on the
# sides, bottom at 1.5 draft)
cv_hx = L_2*abs(np.cos(head_ang)) + 0.5*B*abs(np.sin(head_ang)) + args.cvClear
cv_hy = L_2*abs(np.sin(head_ang)) + 0.5*B*abs(np.cos(head_ang)) + args.cvClear
cv_depth = draft + args.cvClear
hf.console(f"Control volume: x +-{cv_hx:.3f} m, y +-{cv_hy:.3f} m about the hull,"
           f" down to {-cv_depth:.3f} m (clearance {args.cvClear} m)")
update_file('xmin', xbody - cv_hx, path=tsdpath)
update_file('xmax', xbody + cv_hx, path=tsdpath)
update_file('ymin', -cv_hy, path=tsdpath)
update_file('ymax', cv_hy, path=tsdpath)
update_file('zmin', -cv_depth, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

# -- Modify topoSetDict_2
tsdpath = os.path.join('system', 'topoSetDict_2')
hf.console("Modifying topoSetDict_2")
l1 = 3.0
update_file('xmin', xbody - l1*L_2, path=tsdpath)
update_file('xmax', xbody + l1*L_2, path=tsdpath)
update_file('ymin', -l1*L_2, path=tsdpath)
update_file('ymax', l1*L_2, path=tsdpath)
update_file('zmin', -3.0*draft, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)







#### ----- REFINEMENTS ---------- ####
# -- Modify topoSetDict.1
update_file('ybox', ydamp, path='system/topoSetDict.1')
update_file('boxstart', xsponge, path='system/topoSetDict.1')
update_file('boxend', xdamp, path='system/topoSetDict.1')
update_file('zbox', max(3.0*draft, 0.5*lam), path='system/topoSetDict.1')

# -- Modify topoSetDict.2
update_file('ybox', 0.8*ymax, path='system/topoSetDict.2')
update_file('boxstart', 0, path='system/topoSetDict.2')
update_file('boxend', xdamp, path='system/topoSetDict.2')
update_file('zbox', max(2.5*draft, 0.25*lam), path='system/topoSetDict.2')

# -- Modify topoSetDict.3
l3 = max(4.0*L_2, lam)
update_file('ybox', l3, path='system/topoSetDict.3')
update_file('boxstart', xbody - l3, path='system/topoSetDict.3')
update_file('boxend', xbody + l3, path='system/topoSetDict.3')
update_file('zbox', max(2.0*draft, 0.2*lam), path='system/topoSetDict.3')

# -- Modify topoSetDict.4
l4 = 2.5
update_file('ybox', l4*L_2, path='system/topoSetDict.4')
update_file('boxstart', xbody - l4*L_2, path='system/topoSetDict.4')
update_file('boxend', xbody + l4*L_2, path='system/topoSetDict.4')
update_file('zbox', max(1.8*draft, 0.18*lam), path='system/topoSetDict.4')

# -- Modify topoSetDict.5 and .6: like the control volume, rectangles with their
# sides along the mesh axes around the rotated hull, with a clearance of
# ref_clear*L on each side (.6 keeps the control volume 0.2 L inside it)
# -- Modify topoSetDict.5
ref_clear = 0.5
ref_hx = L_2*abs(np.cos(head_ang)) + 0.5*B*abs(np.sin(head_ang)) + ref_clear*2*L_2
ref_hy = L_2*abs(np.sin(head_ang)) + 0.5*B*abs(np.cos(head_ang)) + ref_clear*2*L_2
update_file('ybox', ref_hy, path='system/topoSetDict.5')
update_file('boxstart', xbody - ref_hx, path='system/topoSetDict.5')
update_file('boxend', xbody + ref_hx, path='system/topoSetDict.5')
update_file('zbox', max(1.6*draft, 0.16*lam), path='system/topoSetDict.5')


# -- Modify topoSetDict.6
ref_clear = 0.3
ref_hx = L_2*abs(np.cos(head_ang)) + 0.5*B*abs(np.sin(head_ang)) + ref_clear*2*L_2
ref_hy = L_2*abs(np.sin(head_ang)) + 0.5*B*abs(np.cos(head_ang)) + ref_clear*2*L_2
update_file('ybox', ref_hy, path='system/topoSetDict.6')
update_file('boxstart', xbody - ref_hx, path='system/topoSetDict.6')
update_file('boxend', xbody + ref_hx, path='system/topoSetDict.6')
update_file('zbox', max(1.4*draft, 0.12*lam), path='system/topoSetDict.6')

# -- Modify topoSetDict.7
l4 = 1.4
update_file('ybox', l4*L_2, path='system/topoSetDict.7')
update_file('boxstart', xbody - l4*L_2, path='system/topoSetDict.7')
update_file('boxend', xbody + l4*L_2, path='system/topoSetDict.7')
update_file('zbox', max(1.2*draft, 0.1*lam), path='system/topoSetDict.7')

# Refinement boxes nested in depth: the finest (.6) below the control volume by
# zMargin, each coarser one at least zGap deeper than the next finer one.  Equal
# depths give a two-level jump across one face, which snappyHexMesh (hexRef8,
# 'Dangling coarse cells') cannot handle: the sphere patch is then not created.
zGap = max(0.15, 0.04*lam)
z7 = max(1.2*draft, 0.1*lam)
if lam >= 3.5:   # box .7 (refined from lambda 3.5, helpers.py) encloses the control volume
    z7 = max(z7, cv_depth + args.zMargin)
z6 = max(1.4*draft, 0.12*lam, cv_depth + args.zMargin, z7 + zGap)
z5 = max(1.6*draft, 0.16*lam, z6 + zGap)
z3 = max(2.0*draft, 0.2*lam, z5 + zGap)
z1 = max(3.0*draft, 0.5*lam, z3 + zGap)
for _i, _z in ((1, z1), (3, z3), (5, z5), (6, z6), (7, z7)):
    update_file('zbox', _z, path=f'system/topoSetDict.{_i}')
hf.console(f"Refinement box depths: .1 {z1:.3f}, .3 {z3:.3f}, .5 {z5:.3f}, .6 {z6:.3f}, .7 {z7:.3f} m"
           f" (control volume {cv_depth:.3f} m)")

# For lambda >= 3.5, where helpers.py also refines box .7, the finest box encloses
# the control volume with zMargin to spare on the sides as well, and .6, .5 grow
# around it, each at least zGap wider than the next finer one (sphere_diag test J,
# 2026-10-15: with .7 at +-0.7 m the control volume lay just outside the finest
# level and the lambda 3.9 m surge drift was 4-6 % low).
if lam >= 3.5:
    _b = lambda c: (L_2*abs(np.cos(head_ang)) + 0.5*B*abs(np.sin(head_ang)) + c,
                    L_2*abs(np.sin(head_ang)) + 0.5*B*abs(np.cos(head_ang)) + c)
    h7x, h7y = cv_hx + args.zMargin, cv_hy + args.zMargin
    h6x, h6y = (max(a, b + zGap) for a, b in zip(_b(0.3*2*L_2), (h7x, h7y)))
    h5x, h5y = (max(a, b + zGap) for a, b in zip(_b(0.5*2*L_2), (h6x, h6y)))
    for _i, (_hx, _hy) in ((5, (h5x, h5y)), (6, (h6x, h6y)), (7, (h7x, h7y))):
        update_file('ybox', _hy, path=f'system/topoSetDict.{_i}')
        update_file('boxstart', xbody - _hx, path=f'system/topoSetDict.{_i}')
        update_file('boxend', xbody + _hx, path=f'system/topoSetDict.{_i}')
    hf.console(f"Boxes .7 +-{h7x:.3f}, .6 +-{h6x:.3f}, .5 +-{h5x:.3f} m around the control volume +-{cv_hx:.3f} m")
"""
# -- Modify topoSetDict.6
l2 = 1.15
p1 = np.array([- l2*L_2, l2*B, 0])
p2 = np.array([l2*L_2, l2*B, 0])
p3 = np.array([- l2*L_2, -l2*B, 0])

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
k_z = - max(1.4*depth, 0.05*lam)

update_file('xbox', p1_rot[0], path='system/topoSetDict.6')
update_file('ybox', p1_rot[1], path='system/topoSetDict.6')
update_file('zbox', 0.0, path='system/topoSetDict.6')
update_file('i_x', i_[0], path='system/topoSetDict.6')
update_file('i_y', i_[1], path='system/topoSetDict.6')
update_file('j_x', j_[0], path='system/topoSetDict.6')
update_file('j_y', j_[1], path='system/topoSetDict.6')
"""

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
update_file('writeInterval', f'{100.0:.2f}', path=cdpath)
update_file('cvPoint', f'({xbody} 0 { -draft/2})', path=cdpath)
update_file('CofR', f'({xbody} 0 0)', path=cdpath)
# psiBar (middleFieldForm): free-surface forcing only within 2.5 lambda of the
# hull.  The sponges start 5 lambda from it (inlet xsponge, outlet xdamp, sides
# ydamp), where the damped disturbance would give spurious forcing.  No-op
# unless the controlDict sets psiBarRadius.
update_file('psiBarRadius', f'{2.5*lam:.4g}', path=cdpath)

# -- Modify rigidBodyMotionProperties
rgpath = os.path.join('constant', 'bodyMotionProperties')
update_file('xG', f'{xbody:.4}', path = rgpath)
update_file('Lpp', f'{2*L_2:.6}', path = rgpath)
update_file('beam', f'{B:.6}', path = rgpath)
update_file('heading', f'{np.pi - head_ang:.12}', path = rgpath)

# -- Modify linMotions.py
rgpath = 'linMotions.py'
update_file('lam', f' = {lam:.2}', path = rgpath,  endl='')
update_file('zeta', f' = {zeta0:.4}', path = rgpath,  endl='')
update_file('heading', f' = {heading:.12}', path = rgpath,  endl='')
update_file('Ucur', f' = {Ucur:.4}', path = rgpath,  endl='')

# -- Modify 0.orig/PhiS  (u = +grad(Phi): far field is +U0*x)
phicurpath = os.path.join('0.orig', 'PhiS')
update_file('Ucur', f'{Ucur:.4}', path = phicurpath)
update_file('heading', f'{-heading:.12}', path = phicurpath)
update_file('swaySpeed', f'{np.tan(np.radians(beta)) * Ucur:.12}', path = phicurpath)

# -- translate surface
subprocess.run(['surfaceTransformPoints', '-rotate', f'((-1 0 0) ({-np.cos(head_ang)} {np.sin(head_ang)} 0))', 'constant/triSurface/' + hull + '.stl', 'constant/triSurface/' + hull + '_rotated.stl'])
subprocess.run(['surfaceTransformPoints', '-translate', f'({xbody} 0 0)', 'constant/triSurface/' + hull + '_rotated.stl', 'constant/triSurface/' + hull + '_moved.stl'])

# -- Ready to run
if meshFlag:
    helpers.mesh(lam)
if runFlag:
    rc = helpers.run_case(Nproc)
    print("shipFlow finished rc=", rc)
if postFlag:
    helpers.post_process()

