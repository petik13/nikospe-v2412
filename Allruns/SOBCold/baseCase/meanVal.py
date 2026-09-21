import numpy as np
from scipy.interpolate import UnivariateSpline
import glob
import io
import os


def load(pattern):
    """Load OpenFOAM tabular output, tolerating bracketed vectors."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no files matching {pattern}")
    rows = []
    for f in files:
        raw = open(f, errors="ignore").read().replace("(", " ").replace(")", " ")
        d = np.atleast_2d(np.loadtxt(io.StringIO(raw), comments="#"))
        if d.size:
            rows.append(d)
    if not rows:
        raise ValueError(f"no data in {pattern}")
    d = np.vstack(rows)
    d = d[np.argsort(d[:, 0])]
    _, keep = np.unique(d[:, 0], return_index=True)
    return d[keep]


# -- Load data

heading	 = 30.0
head_ang =  heading * np.pi / 180.0

rho = 1000
zeta	 = 0.025
g = 9.81
lam	 = 2.5
k = 2*np.pi/lam
Ucur	 = 0.3086
omega = np.sqrt(g*k) + Ucur * k * np.cos(head_ang)
scale = 32.0
Lpp = 190/scale
B = 32.2/scale
denom = rho*g*zeta**2*B**2/Lpp
denomMom = denom*Lpp
T = 2*np.pi/omega
print("lam/Lpp = ", lam/Lpp)

fnameForce = f"postProcessing/meanLoads/*/force.dat"
fnameMoment = f"postProcessing/meanLoads/*/moment.dat"

data = load(fnameForce)
dataMom = load(fnameMoment)

# columns: time, total_(x y z), surface_(x y z), elevation_(x y z), strip_(x y z)
time = data[:,0]
Fx = data[:,1]
Fy = data[:,2]
Mz = dataMom[:,3]

i1 = -5
i2 = -1

# -- Find upward crossing indices
upward_crossings = np.where((Fx[:-1] > 0) & (Fx[1:] <= 0))[0]
idx0 = upward_crossings[i1]
idx1 = upward_crossings[i2]
mean_Fx = np.mean(Fx[idx0:idx1])
# upward_crossings = np.where((Fy[:-1] > 0) & (Fy[1:] <= 0))[0]
# idx0 = upward_crossings[i1]
# idx1 = upward_crossings[i2]
mean_Fy = np.mean(Fy[idx0:idx1])
# upward_crossings = np.where((Mz[:-1] > 0) & (Mz[1:] <= 0))[0]
# idx0 = upward_crossings[i1]
# idx1 = upward_crossings[i2]
mean_Mz = np.mean(Mz[idx0:idx1])

meanFx_p = mean_Fx * np.cos(head_ang) - mean_Fy * np.sin(head_ang)
meanFy_p = mean_Fx * np.sin(head_ang) + mean_Fy * np.cos(head_ang)

expFolder = os.path.expanduser('~/expWadRes/')
wamDat = np.loadtxt(expFolder + 'wadDat30.dat')
expFx = np.loadtxt(expFolder + 'expFx30U0.dat')
expFy = np.loadtxt(expFolder + 'expFy30U0.dat')
expMz = np.loadtxt(expFolder + 'expMz30U0.dat')
expFxU = np.loadtxt(expFolder + 'expFx30U0387.dat')
expFyU = np.loadtxt(expFolder + 'expFy30U0387.dat')
expMzU = np.loadtxt(expFolder + 'expMz30U0387.dat')


meanFx_p_ND = meanFx_p/denom
meanFy_p_ND = meanFy_p/denom
meanMz_ND = mean_Mz/denomMom

with open('../sphereFixed_U' + str(Ucur) + '.dat', 'a') as f:
    f.write(f"{lam/Lpp} {omega**2*B/2/g} {meanFx_p_ND} {meanFy_p_ND} {meanMz_ND} \n")

print("omega**2*R/g = ", omega**2*B/2/g)
print("meanFx_p_ND = ", meanFx_p_ND)
print("meanFy_p_ND = ", meanFy_p_ND)
print("meanMz_ND = ", meanMz_ND)

print('tau = ', omega*Ucur/9.81)
