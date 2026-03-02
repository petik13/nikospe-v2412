from matplotlib.cbook import ls_mapper
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

fname = 'postProcessing/patchProbes1/0/zeta'

with open(fname) as f:
    cleaned = (line.replace("(", "").replace(")", "") for line in f)
    data = np.loadtxt(cleaned)

time = data[:,0]
z1   = data[:,3]

fig, ax = plt.subplots()
ax.plot(time, z1, label='z1')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Surface Elevation (m)')

# find peaks
peaks, _ = find_peaks(z1)

# extract peak values
peak_values = z1[peaks]

# take last 10 peaks
last10 = peak_values[-10:]

# compute statistics
mean_val = np.mean(last10)
std_val  = np.std(last10)

print("Last 10 peak values:", last10)
print("Mean:", mean_val)
print("Standard deviation:", std_val)

heading = 30
head_ang = np.radians(30)
lam = 2.0
k = 2 * np.pi / lam
Ucur = 0.387
omega = np.sqrt(9.81 * k ) + Ucur * k * np.cos(head_ang)
T = 2 * np.pi / omega
steepness = 1/60
zeta = steepness * lam * 0.5
print("Calculated period from wavelength and current:", T, "seconds")
Tnum = time[peaks][-1] - time[peaks][-2]
print("Calculated period from last two peaks:", Tnum, "seconds")

# optional: visualize
ax.plot(time, z1)
ax.plot(time[peaks], z1[peaks], "ro")

# plot horizontal line at zeta = 0.0166
ax.axhline(y=zeta, color='g', linestyle='--', label=f'zeta = {zeta:.4f} m')

plt.show()
