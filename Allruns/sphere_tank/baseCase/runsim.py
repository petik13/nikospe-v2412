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
# Hemisphere in the NTH towing tank of Zhao, Faltinsen, Krokstad & Aanesland
# (1988), "Wave-current interaction effects on large-volume structures", BOSS'88:
# tank breadth 10.5 m, depth 10 m (length 85 m), hemisphere D = 1 m on the tank
# centreline, restrained, wave height H = 0.1 m, current simulated by towing the
# model against the waves (U > 0: current along the waves, omega_e > omega_0).
#
# Modelled: the side walls (front/back: walls, no damping zone, zero normal
# gradient for PhiD, PhiS and psibar) and the depth.  Not modelled: the tank
# length, wave maker and beach -- as in Zhao et al.'s own tank-wall theory
# (13 image singularities on each side, infinitely long tank), the run is in
# the body frame with damping zones at the inlet and outlet.  The side walls
# are only consistent with waves along the tank: heading 0 (head sea, as in
# the paper) or 180.
#
# Give the wave either as --lam [m] or, as in the paper, as the encounter
# period --Te [s] (lambda then follows from omega_e = omega_0 + k U, finite depth).
parser = argparse.ArgumentParser()
parser.add_argument("--lam", type=float, default=None, help="wave length [m]")
parser.add_argument("--Te", type=float, default=None,
                    help="encounter period [s] (Zhao et al. Figs 7-9: 1.14 ... 1.57 s)")
parser.add_argument("--heading", type=float, default=0.0)
parser.add_argument("--Ucur", type=float, required=True,
                    help="towing speed [m/s] (Zhao et al.: -0.1, 0, 0.1, 0.2)")
parser.add_argument("--H", type=float, default=0.1, help="wave height [m] (Zhao et al.: 0.1)")
# Control-volume placement (sphere_diag tests A-G, 2026-10-15): the midfield drift
# force is converged at the base resolution only on a control surface ~8 cells or
# more from the hull and entirely inside the finest refinement box, also in depth.
# 0.1 L from the hull (3 cells) over-predicted the surge drift by 8 %.
parser.add_argument("--cvClear", type=float, default=0.25,
                    help="clearance of the control volume around the hull, sides and below the keel [m]")
parser.add_argument("--zMargin", type=float, default=0.1,
                    help="finest refinement boxes at least this much deeper than the control volume [m]")
args = parser.parse_args()

TANK_B = 10.5        # tank breadth [m]
TANK_DEPTH = 10.0    # water depth [m]
TANK_L = 85.0        # tank length [m] (not modelled, see above)

if min(abs(args.heading % 360.0), abs(args.heading % 360.0 - 180.0)) > 1e-9:
    raise SystemExit("sphere_tank: the side walls need waves along the tank, heading 0 or 180")
if (args.lam is None) == (args.Te is None):
    raise SystemExit("sphere_tank: give exactly one of --lam and --Te")


def lam_from_Te(Te, U, heading_deg, h=TANK_DEPTH, g=9.81):
    """Wave length whose encounter period is Te: omega_e = sqrt(g k tanh(k h)) + k U cos(heading)."""
    c = np.cos(np.radians(heading_deg))
    f = lambda k: np.sqrt(g*k*np.tanh(k*h)) + k*U*c - 2*np.pi/Te
    lo, hi = 1e-4, 1e3
    if U*c < 0:
        # stream against the waves: omega_e(k) peaks at k = g/(4 U^2) (deep
        # water); the physical root is the one below the peak
        hi = min(hi, g/(4*(U*c)**2))
    if f(lo)*f(hi) > 0:
        raise SystemExit(f"sphere_tank: no wave with encounter period {Te} s at U = {U} m/s")
    for _ in range(200):                                   # bisection
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if f(lo)*f(mid) > 0 else (lo, mid)
    return 2*np.pi/(0.5*(lo + hi))


# -- Define parameters
hull = 'sphere'
heading = args.heading
lam = args.lam if args.lam is not None else lam_from_Te(args.Te, args.Ucur, args.heading)
steepness = args.H/lam
zeta0 = 0.5*args.H
hf.console(f"Tank case: lambda {lam:.5g} m, H {args.H} m (H/lambda {steepness:.4g}),"
           f" U {args.Ucur} m/s, heading {heading} deg")
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

# Side walls at +-B/2 with no damping zone (ySide far outside the domain)
ydamp = 1000.0
Lydamp = 2 * lam
ymin = -0.5*TANK_B
ymax = 0.5*TANK_B
zmin = -TANK_DEPTH
zmax = 0.0
# depth of the open-water sphere case, for its vertical mesh spacing near the surface
zmin_ref = min(-0.5*lam, -7*draft)

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
# Vertical: the open-water case's cells near the surface (Nz_ref cells over
# -zmin_ref, top/bottom ratio 0.2), continued downwards with the same expansion
# to the tank floor
Nz_ref = int(-zmin_ref / lam_mesh * discZ)
r_ref = 0.2**(1.0/max(Nz_ref - 1, 1))                   # cell-to-cell ratio, upwards
dz_top = -zmin_ref*(1 - 1/r_ref)/(1 - (1/r_ref)**Nz_ref)  # top cell of the reference
Nz, depth_sum = 0, 0.0
while depth_sum < -zmin:
    depth_sum += dz_top/r_ref**Nz
    Nz += 1
zGrading = r_ref**(Nz - 1)
hf.console(f"Vertical mesh: {Nz} cells to {zmin} m, top cell {dz_top:.4g} m, grading {zGrading:.4g}")

# Timestep calculation
k = 2 * np.pi / lam
g = 9.81
omega = np.sqrt(g * k)
period = 2 * np.pi / omega
celerity = omega / k
Cgroup = 0.5*celerity
T = 2*np.pi/omega
endTime = rampperiod*T + 0.5*(xbody + xbody)/Cgroup
# the waves the walls reflect need time to cross the tank and back a few times
endTime += 3*TANK_B/Cgroup
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
update_file('zGrading', f'{zGrading:.6g}')
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
# -- Modify topoSetDict.1 (the whole tank breadth: the walls reflect, so the
# waves matter everywhere between them)
update_file('ybox', 0.5*TANK_B, path='system/topoSetDict.1')
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
# psiBar (middleFieldForm): free-surface forcing only within psiBarRadius of the
# hull, at least 2.5 lambda and at least out to the side walls (no side sponge
# here); the inlet and outlet sponges start 5 lambda from it, where the damped
# disturbance would give spurious forcing, so it stays below 4 lambda.
psi_radius = min(max(2.5*lam, 0.5*TANK_B), 4.0*lam)
if psi_radius < 0.5*TANK_B:
    hf.console(f"WARNING: psiBarRadius {psi_radius:.3g} m does not reach the walls"
               f" ({0.5*TANK_B} m): lambda {lam:.3g} m is short for this tank layout")
update_file('psiBarRadius', f'{psi_radius:.4g}', path=cdpath)

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

