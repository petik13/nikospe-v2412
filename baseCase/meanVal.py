import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import os
# -- Load data

heading	 = 60.0
head_ang =  heading * np.pi / 180.0

rho = 1000
zeta	 = 0.05
g = 9.81
lam	 = 6.0
k = 2*np.pi/lam
Ucur	 = 0.0
omega = np.sqrt(g*k) + Ucur * k * np.cos(head_ang)
denom = zeta**2
T = 2*np.pi/omega


fnameForce = f"postProcessing/MWLforces/0/force.dat"
fnameMoment = f"postProcessing/MWLforces/0/moment.dat"

data = np.loadtxt(fnameForce, skiprows=1)
dataMom = np.loadtxt(fnameMoment, skiprows=1)

time = data[:,0]
Fx = data[:,1]
Fy = data[:,2]
Mz = dataMom[:,3]

i1 = -10
i2 = -1
# -- Find upward crossing indices
upward_crossings = np.where((Fx[:-1] > 0) & (Fx[1:] <= 0))[0]
idx0 = upward_crossings[i1]
idx1 = upward_crossings[i2]
mean_Fx = np.mean(Fx[idx0:idx1])
print("mean_Fx:", mean_Fx)
#upward_crossings = np.where((Fy[:-1] > 0) & (Fy[1:] <= 0))[0]
#idx0 = upward_crossings[i1]
#idx1 = upward_crossings[i2]
mean_Fy = np.mean(Fy[idx0:idx1])
print("mean_Fy:", mean_Fy)
#upward_crossings = np.where((Mz[:-1] > 0) & (Mz[1:] <= 0))[0]
#idx0 = upward_crossings[i1]
#idx1 = upward_crossings[i2]
mean_Mz = np.mean(Mz[idx0:idx1])
print("mean_Mz", mean_Mz)

meanFx_p = mean_Fx * np.cos(head_ang) - mean_Fy * np.sin(head_ang)
meanFy_p = mean_Fx * np.sin(head_ang) + mean_Fy * np.cos(head_ang)

expFolder = os.path.expanduser('~/expWadRes/')
wamDat = np.loadtxt(expFolder + 'wadDat60.dat')
expFx = np.loadtxt(expFolder + 'expFx60U0.dat')
expFy = np.loadtxt(expFolder + 'expFy60U0.dat')
expMz = np.loadtxt(expFolder + 'expMz60U0.dat')
expFxU = np.loadtxt(expFolder + 'expFx60U0387.dat')
expFyU = np.loadtxt(expFolder + 'expFy60U0387.dat')
expMzU = np.loadtxt(expFolder + 'expMz60U0387.dat')


meanFx_p_ND = meanFx_p/denom
meanFy_p_ND = meanFy_p/denom
meanMz_ND = mean_Mz/denom

with open('../res30U.dat', 'a') as f:
    f.write(f"{T} {meanFx_p_ND} {meanFy_p_ND} {meanMz_ND}\n")

figs, axs = plt.subplots(1, 3, figsize=(18, 5))
axs[0].plot(wamDat[:,0], -wamDat[:,1], 'k-', label='Wadam')
axs[0].plot(T, -meanFx_p_ND, 'rx', label='Simulation Result')
axs[0].plot(expFx[:,0], expFx[:,1], 'ro', label='Experiment')
axs[0].plot(expFxU[:,0], expFxU[:,1], 'bo', label='Experiment U=0.387')
axs[0].set_xlabel('T [s]')
axs[0].set_ylabel('$F_x/ζ²$')
axs[0].legend()
print(meanFx_p_ND)

axs[1].plot(wamDat[:,0], -wamDat[:,2], 'k-', label='Wadam')
axs[1].plot(T, meanFy_p_ND, 'rx', label='Simulation Result')
axs[1].plot(expFy[:,0], expFy[:,1], 'ro', label='Experiment')
axs[1].plot(expFyU[:,0], expFyU[:,1], 'bo', label='Experiment U=0.387')
axs[1].set_xlabel('T [s]')
axs[1].set_ylabel('$F_y/ζ²$')
axs[1].legend()
print(meanFy_p_ND)

axs[2].plot(wamDat[:,0], wamDat[:,3], 'k-', label='Wadam')
axs[2].plot(T, meanMz_ND, 'rx', label='Simulation Result')
axs[2].plot(expMz[:,0], -expMz[:,1], 'ro', label='Experiment')
axs[2].plot(expMzU[:,0], -expMzU[:,1], 'bo', label='Experiment U=0.387')
axs[2].set_xlabel('T [s]')
axs[2].set_ylabel('$M_z/ζ²$')
axs[2].legend()
print(meanMz_ND)

print('tau = ', omega*Ucur/9.81)
plt.show()
