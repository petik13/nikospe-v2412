import numpy as np
import matplotlib.pyplot as plt

# -- Load data
fname = "postProcessing/MWLforces/0/force.dat"
data = np.loadtxt(fname, skiprows=1)
time = data[:,0]
Fx = data[:,1]

# -- Find upward crossing indices
upward_crossings = np.where((Fx[:-1] > 0) & (Fx[1:] <= 0))[0]
idx0 = upward_crossings[-2]
idx1 = upward_crossings[-1]

idx0 = np.where(time >= 50.0)[0][0]
idx1 = np.where(time >= 59.9)[0][0]

# -- Plotting
fig, ax = plt.subplots()
mean_Fx = np.mean(Fx[idx0:idx1])
ax.hlines(mean_Fx, time[idx0], time[idx1], colors='r', linestyles='dashed', label=f'Mean Fx after t= {time[idx0]:.2f}s')
ax.plot(time[idx0:idx1], Fx[idx0:idx1], 'b-', label=f'Fx after t= {time[idx0]:.2f}s')
print(f"Mean Fx after t={time[idx0]:.2f}s: {mean_Fx:.2f} N")
ax.set_xlabel('Time')
ax.set_ylabel('Force')
ax.legend()


expDat = np.loadtxt('expDat.dat')
rho = 1000
zeta	 = 0.0207
D = 1
R = D/2
g = 9.81
lam	 = 6.0
k = 2*np.pi/lam
omega = np.sqrt(g*k)

omND = omega**2*R/g
fND = -mean_Fx/(0.5*rho*g*zeta**2*D)
fig, ax = plt.subplots()
ax.plot(expDat[:,0], expDat[:,1], 'o', label='Experimental Data')
ax.plot(omND, fND, 'rx', label='Simulation Result')
ax.set_xlabel('ω²R/g')
ax.set_ylabel('F/(0.5ρgζ²D)')
ax.legend()
plt.show()
