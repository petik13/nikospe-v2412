"""
Circular motion test (CMT) in regular waves with manFlowPrescribed.

The hull stays aligned with the mesh (bow at -x, heading 0 geometry); the ship
motion (u, v, r) and the rotating incident wave are imposed through
constant/prescribedMotion.  Meshing is the same as runsim.py of the SOBC
sweep.

usage (from a copy of this baseCase):
    python3 runCMT.py --lam 3.2 --u 0.318 --v -0.097 --r 0.08
    python3 runCMT.py --lam 3.2 --u 0.33 --v 0 --r 0        # regression vs manFlow
options:
    --psi0 0            initial heading [deg]
    --waveDir           wave propagation direction [deg], default psi0 + 180 (head seas at t = 0)
    --yawOnsetPeriods 10  straight run before the yaw rate starts, in encounter
                        periods at the initial heading (the motions settle)
    --yawRampPeriods 2  duration of the (half-cosine) yaw-rate ramp, same units
    --endTime           default: r != 0: one full turn of mean loads from the onset,
                        onset + ramp/2 + 2 pi/|r|, plus the half-width of the
                        mean-load filter in cmtPost.py (2 Te_max) and a margin
                        of Te_max, so that results.csv reaches the end-of-turn row;
                        r = 0: wave ramp + 0.5*(2 xbody)/c_g, as runsim.py
    --nproc 52 --procD 8 6 1
    --noMesh --noRun --noPost

u, v, r are MMG quantities at midship (x forward, y starboard, r > 0 to
starboard).  The r = 0 run is equivalent to the manFlow run with
--Ucur u and sway speed -v (see constant/prescribedMotion).

Mean loads: meanLoadsRot (control volume +-1.2 L/2) and meanLoadsRot2
(+-2 L/2, topoSetDict_2) are the rotating-frame midfield; they must agree.
"""

import glob
import numpy as np
import shutil
import subprocess
import helperFuncs as hf
import os
import argparse
import helpers

parser = argparse.ArgumentParser()
parser.add_argument("--lam", type=float, required=True, help="wavelength [m]")
parser.add_argument("--u", type=float, required=True, help="surge [m/s]")
parser.add_argument("--v", type=float, default=0.0, help="sway at midship [m/s], starboard +")
parser.add_argument("--r", type=float, default=0.0, help="yaw rate [rad/s], starboard +")
parser.add_argument("--psi0", type=float, default=0.0, help="initial heading [deg]")
parser.add_argument("--waveDir", type=float, default=None, help="wave propagation direction [deg]")
parser.add_argument("--yawOnsetPeriods", type=float, default=10.0,
                    help="straight run before the yaw-rate ramp [encounter periods]")
parser.add_argument("--yawRampPeriods", type=float, default=2.0,
                    help="duration of the yaw-rate ramp [encounter periods]")
parser.add_argument("--endTime", type=float, default=None)
parser.add_argument("--nproc", type=int, default=52)
parser.add_argument("--procD", type=int, nargs=3, default=[8, 6, 1])
parser.add_argument("--noMesh", action="store_true")
parser.add_argument("--noRun", action="store_true")
parser.add_argument("--noPost", action="store_true")
args = parser.parse_args()

logger = hf.setup_logger('log')
hf.console("Starting CMT workflow")
subprocess.run(['./clean.sh'])
if os.path.isdir('0'):
    shutil.rmtree('0')
shutil.copytree('0.orig', '0')

# -- Define parameters (as runsim.py)
hull = 'SOBC'
steepness = 1/50
lam = args.lam
zeta0 = 0.5*steepness*lam
scale = 32.0
L_2 = 190.0/scale/2
B = 32.2/scale
draft = 11.0/scale
rampperiod = 3.0
Co = 0.2
Nproc = args.nproc
procD = args.procD

u, v, r = args.u, args.v, args.r
psi0 = args.psi0
waveDir = psi0 + 180.0 if args.waveDir is None else args.waveDir

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
Nref = 2
if lam >= 1.5:
    Nref = 3
if lam >= 3.0:
    Nref = 4
if lam >= 4.5:
    Nref = 5

lam_mesh = lam
Nx = int(xmax / lam_mesh * discX)
Ny = int((ymax - ymin) / lam_mesh * discY)
Nz = int(-zmin / lam_mesh * discZ)

# -- Frame velocity V_S = (-u, v) + (r dy, -r dx) relative to midship
def VS_mag(x, y):
    dx, dy = x - xbody, y
    return np.hypot(-u + r*dy, v - r*dx)

l3 = max(2.2*L_2, lam)                          # finest refinement box, half size
V_near = max(VS_mag(xbody + sx*l3, sy*l3) for sx in (-1, 1) for sy in (-1, 1))
V_far = max(VS_mag(x, y) for x in (xmin, xmax) for y in (ymin, ymax))

# -- Timestep: as runsim.py with the ship speed replaced by the largest frame
#    velocity, checked separately on the fine (near) and coarse (far) cells
k = 2 * np.pi / lam
g = 9.81
omega = np.sqrt(g * k)
celerity = omega / k
Cgroup = 0.5*celerity
T = 2*np.pi/omega
dt_near = Co * lam_mesh / ((celerity + V_near) * discX)/(Nref + 1)
dt_far = Co * lam_mesh / ((celerity + V_far) * discX)
deltaT = min(dt_near, dt_far)

# -- Yaw-rate onset, in encounter periods at the initial heading (finite-depth
#    dispersion, as the solver).  Wave direction in mesh axes e = (-cos a, sin a),
#    a = chi - psi0; frame velocity V0 = (-u, v); omega_e = |omega - k e.V0|.
omega_h = np.sqrt(g*k*np.tanh(k*(-zmin)))
a0 = np.radians(waveDir - psi0)
eV0 = np.cos(a0)*u + np.sin(a0)*v
Te0 = 2*np.pi/abs(-omega_h + k*eV0)
tOn = args.yawOnsetPeriods*Te0
tRamp = args.yawRampPeriods*Te0

# Longest encounter period over a full turn.  cmtPost.py smooths the mean loads
# with a centred Gaussian, sigma = 0.5 Te_max, which needs 4 sigma = 2 Te_max of
# data on either side: before the onset (the straight run) and after the turn.
al = np.linspace(0.0, 2*np.pi, 721)
we_turn = np.abs(-omega_h + k*(np.cos(al)*u + np.sin(al)*v))
Te_max = 2*np.pi/max(we_turn.min(), 0.05*omega_h)
tFilter = 2.0*Te_max

if args.endTime is not None:
    endTime = args.endTime
elif r != 0:
    # After the half-cosine ramp the heading is psi0 + r (t - tOn - tRamp/2):
    # one full turn of mean loads from the onset, plus the filter half-width
    endTime = tOn + 0.5*tRamp + 2*np.pi/abs(r) + tFilter + Te_max   # + Te_max margin: the end row
else:
    endTime = rampperiod*T + 0.5*(xbody + xbody)/Cgroup

hf.console(f"CMT: u {u} m/s  v {v} m/s  r {r} rad/s  (beta {np.degrees(np.arctan2(-v, u)):.1f} deg,"
           f" turn radius {np.hypot(u, v)/r if r else np.inf:.3g} m)")
hf.console(f"     encounter period at the start {Te0:.3f} s; yaw onset at {tOn:.2f} s"
           f" ({args.yawOnsetPeriods:g} periods), ramp {tRamp:.2f} s ({args.yawRampPeriods:g} periods)")
if r != 0:
    hf.console(f"     heading change over the run "
               f"{np.degrees(r*max(endTime - tOn - 0.5*tRamp, 0)):.1f} deg;"
               f" mean loads from the onset (heading change 0) to "
               f"{endTime - tFilter:.2f} s (filter half-width {tFilter:.2f} s)")
    if tOn < tFilter + rampperiod*T:
        hf.console(f"     WARNING: the straight run before the onset ({tOn:.2f} s) is shorter"
                   f" than the wave ramp + filter half-width ({tFilter + rampperiod*T:.2f} s);"
                   f" the mean loads at the onset will be affected by the start-up")
hf.console(f"     max |V_S|: near {V_near:.3f} m/s, far {V_far:.3f} m/s")
hf.console(f"     deltaT = {deltaT:.6f} s (near {dt_near:.6f}, far {dt_far:.6f}),"
           f" endTime = {endTime:.2f} s")

def update_file(name, value, path='system/blockMeshDict', endl=';'):
    if not os.path.isfile(path):
        logger.warning(f"file not found at {path}, skipping update.")
        return

    with open(path, 'r') as f:
        lines = f.readlines()

    def repl(line: str, key: str, value) -> str:
        stripped = line.lstrip()
        if not stripped or stripped.startswith('//') or stripped.startswith('/*'):
            return line
        first = stripped.split(None, 1)[0]
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
update_file('bodyXpos', xbody, path='system/snappyHexMeshDict')

# -- waveConditions: wave and sponges.  currentSpeed/headingAngle/swaySpeed
#    are not used by the rotating-frame solver; set for the record.
wCpath = os.path.join('constant', 'waveConditions')
hf.console("Modifying waveConditions")
update_file('headingAngle', 0.0, path=wCpath)
update_file('steepness', steepness, path=wCpath)
update_file('waveLength', lam, path=wCpath)
update_file('currentSpeed', u, path=wCpath)
update_file('swaySpeed', -v, path=wCpath)
update_file('waterDepth', -zmin, path=wCpath)
update_file('xOutlet', xdamp, path=wCpath)
update_file('LOutlet', Lxdamp, path=wCpath)
update_file('ySide', ydamp, path=wCpath)
update_file('LSide', Lydamp, path=wCpath)
update_file('xInlet', xsponge, path=wCpath)
update_file('LInlet', Lsponge, path=wCpath)
update_file('rampPeriods', rampperiod, path=wCpath)

# -- prescribedMotion
pmpath = os.path.join('constant', 'prescribedMotion')
hf.console("Modifying prescribedMotion")
update_file('u', f'{u:.6g}', path=pmpath)
update_file('v', f'{v:.6g}', path=pmpath)
update_file('r', f'{r:.6g}', path=pmpath)
update_file('psi0', f'{psi0:.6g}', path=pmpath)
update_file('waveDirection', f'{waveDir:.6g}', path=pmpath)
update_file('rotationCentre', f'({xbody} 0 0)', path=pmpath)
update_file('yawOnsetTime', f'{tOn:.6g}', path=pmpath)
update_file('yawRampTime', f'{tRamp:.6g}', path=pmpath)

####### CONTROL VOLUMES: controlZone (+-1.2 L/2), controlZone2 (+-2 L/2) ########
# controlZone2 is kept close to the hull, inside the refined region (refinement
# box 4 is +-2.5 L/2): further out the cells are coarser and the numerical
# errors of the wave field larger, which would spoil the comparison.
tsdpath = os.path.join('system', 'topoSetDict')
l1 = 1.2
update_file('xmin', xbody - l1*L_2, path=tsdpath)
update_file('xmax', xbody + l1*L_2, path=tsdpath)
update_file('ymin', -l1*L_2, path=tsdpath)
update_file('ymax', l1*L_2, path=tsdpath)
update_file('zmin', -1.5*draft, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

tsdpath = os.path.join('system', 'topoSetDict_2')
l1 = 2.0
update_file('xmin', xbody - l1*L_2, path=tsdpath)
update_file('xmax', xbody + l1*L_2, path=tsdpath)
update_file('ymin', -l1*L_2, path=tsdpath)
update_file('ymax', l1*L_2, path=tsdpath)
update_file('zmin', -2.0*draft, path=tsdpath)
update_file('zmax', zmax, path=tsdpath)

#### ----- REFINEMENTS (as runsim.py) ---------- ####
update_file('ybox', ydamp, path='system/topoSetDict.1')
update_file('boxstart', xsponge, path='system/topoSetDict.1')
update_file('boxend', xdamp, path='system/topoSetDict.1')
update_file('zbox', max(3.0*draft, 0.5*lam), path='system/topoSetDict.1')

update_file('ybox', 0.8*ymax, path='system/topoSetDict.2')
update_file('boxstart', 0, path='system/topoSetDict.2')
update_file('boxend', xdamp, path='system/topoSetDict.2')
update_file('zbox', max(2.5*draft, 0.25*lam), path='system/topoSetDict.2')

update_file('ybox', l3, path='system/topoSetDict.3')
update_file('boxstart', xbody - l3, path='system/topoSetDict.3')
update_file('boxend', xbody + l3, path='system/topoSetDict.3')
update_file('zbox', max(2.0*draft, 0.2*lam), path='system/topoSetDict.3')

for i, l4, zf, zl in ((4, 2.5, 1.8, 0.18), (5, 1.4, 1.6, 0.16), (6, 1.3, 1.4, 0.12), (7, 1.2, 1.2, 0.1)):
    p = f'system/topoSetDict.{i}'
    update_file('ybox', l4*L_2, path=p)
    update_file('boxstart', xbody - l4*L_2, path=p)
    update_file('boxend', xbody + l4*L_2, path=p)
    update_file('zbox', max(zf*draft, zl*lam), path=p)

# -- decomposeParDict
dpdpath = os.path.join('system', 'decomposeParDict')
update_file('numberOfSubdomains', Nproc, path=dpdpath)
update_file('n', f'\t({procD[0]} {procD[1]} {procD[2]}) ', path=dpdpath)

# -- controlDict
cdpath = os.path.join('system', 'controlDict')
hf.console("Modifying controlDict")
update_file('deltaT', f'{deltaT:.5f}', path=cdpath)
update_file('endTime', f'{endTime:.2f}', path=cdpath)
update_file('writeInterval', f'{100.0:.2f}', path=cdpath)
update_file('cvPoint', f'({xbody} 0 { -draft/2})', path=cdpath)
update_file('CofR', f'({xbody} 0 0)', path=cdpath)

# -- bodyMotionProperties: hull aligned with the mesh, bow at -x
rgpath = os.path.join('constant', 'bodyMotionProperties')
update_file('xG', f'{xbody:.4}', path=rgpath)
update_file('Lpp', f'{2*L_2:.6}', path=rgpath)
update_file('beam', f'{B:.6}', path=rgpath)
update_file('heading', f'{np.pi:.12}', path=rgpath)

# -- hull surface: heading 0 geometry, translated to xbody
subprocess.run(['surfaceTransformPoints', '-rotate', '((-1 0 0) (-1 0 0))',
                'constant/triSurface/' + hull + '.stl', 'constant/triSurface/' + hull + '_rotated.stl'])
subprocess.run(['surfaceTransformPoints', '-translate', f'({xbody} 0 0)',
                'constant/triSurface/' + hull + '_rotated.stl', 'constant/triSurface/' + hull + '_moved.stl'])


def run_case(Nproc):
    subprocess.run(['rm', '-r', '0'])
    subprocess.run(['cp', '-r', '0.orig', '0'])
    subprocess.run(['topoSet', '-dict', 'system/topoSetDict'])
    subprocess.run(['topoSet', '-dict', 'system/topoSetDict_2'])
    subprocess.run(['renumberMesh', '-overwrite'])
    subprocess.run(['decomposePar'])
    subprocess.run(['foamJob', '-s', '-p', 'renumberMesh', '-overwrite'])

    with open("log", "w") as log:
        proc = subprocess.Popen(
            ['mpirun', '-np', str(Nproc),
             'manFlowPrescribed', '-parallel', '-withFunctionObjects'],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        rc = proc.wait()
    return rc


def clean_processors():
    """Reconstruct the written fields (if any) and remove the processor
    directories, as prepPost.sh does for the straight-course runs.  The
    post-processing only needs postProcessing/."""
    subprocess.run(['reconstructPar'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for d in glob.glob('processor*'):
        shutil.rmtree(d, ignore_errors=True)


def plot_series():
    """Time series of the rotating-midfield loads and the motions
    (timeSeries.png).  Without DISPLAY plotSeries.py writes the file instead
    of opening a window, so the workflow never blocks."""
    env = dict(os.environ, DISPLAY="")
    subprocess.run(['python3', 'plotSeries.py', '--save', '--fo', 'meanLoadsRot'], env=env)


# -- Ready to run
if not args.noMesh:
    helpers.mesh(lam)
if not args.noRun:
    rc = run_case(Nproc)
    print("manFlowPrescribed finished rc=", rc)
    if rc == 0:
        clean_processors()
    else:
        print("run failed: processor directories kept (restart / debugging)")
if not args.noPost:
    subprocess.run(['python3', 'cmtPost.py'])
    plot_series()
