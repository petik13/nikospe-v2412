from math import e
import numpy as np
import matplotlib.pyplot as plt
import os

Ucur	 = 0.0
heading	 = 60.0
head_ang = np.deg2rad(heading)
lam	 = 6.0
k = 2*np.pi/lam
omega = np.sqrt(9.81 * k) + Ucur*k*np.cos(head_ang)
T = 2*np.pi/omega
zeta	 = 0.05

expFolder = os.path.expanduser('~/expWadRes/')

fname = f'postProcessing/rigidBodyMotion/motionTurgut.dat'
data = np.loadtxt(fname, skiprows=1)

time = data[:,0]

idx1 = np.where((np.diff(data[:-1, 1]) > 0) & (np.diff(data[1:, 1]) < 0))[0] + 1
eta1 = np.mean(data[idx1[-6:-1],1])/zeta
idx2 = np.where((np.diff(data[:-1, 2]) > 0) & (np.diff(data[1:, 2]) < 0))[0] + 1
eta2 = np.mean(data[idx2[-6:-1],2])/zeta
idx3 = np.where((np.diff(data[:-1, 3]) > 0) & (np.diff(data[1:, 3]) < 0))[0] + 1
eta3 = np.mean(data[idx3[-6:-1],3])/zeta
idx4 = np.where((np.diff(data[:-1, 4]) > 0) & (np.diff(data[1:, 4]) < 0))[0] + 1
eta4 = np.mean(data[idx4[-6:-1],4])/zeta
idx5 = np.where((np.diff(data[:-1, 5]) > 0) & (np.diff(data[1:, 5]) < 0))[0] + 1
eta5 = np.mean(data[idx5[-6:-1],5])/zeta
idx6 = np.where((np.diff(data[:-1, 6]) > 0) & (np.diff(data[1:, 6]) < 0))[0] + 1
eta6 = np.mean(data[idx6[-6:-1],6])/zeta
# print(eta1)
wadLin = np.loadtxt(expFolder + 'wadLinMotions.dat')

figs, axs = plt.subplots(2, 3, figsize=(18, 10))
axs[0,0].plot(wadLin[:,0], wadLin[:,1], 'k-', label = 'Wadam')
axs[0, 1].plot(wadLin[:,0], wadLin[:,2], 'k-', label = 'Wadam')
axs[0, 2].plot(wadLin[:,0], wadLin[:,3], 'k-', label = 'Wadam')
axs[1,0].plot(wadLin[:,0], wadLin[:,4]*63.65, 'k-', label = 'Wadam')
axs[1,1].plot(wadLin[:,0], wadLin[:,5]*63.65, 'k-', label = 'Wadam')
axs[1,2].plot(wadLin[:,0], wadLin[:,6]*63.65, 'k-', label = 'Wadam')


# plot numerical results
axs[0,0].plot(T, eta1, 'ro', label='Simulation Result')
axs[0,1].plot(T, eta2, 'ro', label='Simulation Result')
axs[0,2].plot(T, eta3, 'ro', label='Simulation Result')
axs[1,0].plot(T, eta4, 'ro', label='Simulation Result')
axs[1,1].plot(T, eta5, 'ro', label='Simulation Result')
axs[1,2].plot(T, eta6, 'ro', label='Simulation Result')

for ax in axs.flat:
    ax.set_xlabel('T [s]')
    ax.legend()
    ax.grid()
axs[0,0].set_ylabel('$η_1/ζ$')
axs[0,1].set_ylabel('$η_2/ζ$')
axs[0,2].set_ylabel('$η_3/ζ$')
axs[1,0].set_ylabel('$η_4$')
axs[1,1].set_ylabel('$η_5$')
axs[1,2].set_ylabel('$η_6$')


# ================= FFT ANALYSIS OF HEAVE =================

# Extract time and heave
t = data[:,0]
heave = data[:,4]

# Remove mean (important!)
heave = heave - np.mean(heave)

# Optionally use only last portion (remove transient)
frac = 0.7   # use last 50%
start = int(len(t)*(1-frac))
t_fft = t[start:]
heave_fft = heave[start:]

# Sampling info
dt = np.mean(np.diff(t_fft))
fs = 1/dt
N = len(heave_fft)

# FFT
fft_vals = np.fft.rfft(heave_fft)
freqs = np.fft.rfftfreq(N, dt)

# Amplitude spectrum (properly scaled)
amp = 2*np.abs(fft_vals)/N

# Ignore zero frequency
freqs = freqs[1:]
amp = amp[1:]

# Convert to period
periods = 1/freqs

# ================= FIND DOMINANT PEAKS =================

# Sort by amplitude
idx_sorted = np.argsort(amp)[::-1]

nPeaks = 5
print("\nDominant frequency peaks:")
print("Rank | Frequency [Hz] | Period [s] | Amplitude")

for i in range(nPeaks):
    idx = idx_sorted[i]
    print(f"{i+1:4d} | {freqs[idx]:.5f} | {periods[idx]:.5f} | {amp[idx]:.5e}")

print("\nEncounter period (theoretical):", T)

# ================= PLOTTING =================

fig, axs = plt.subplots(3,1, figsize=(8,10))

# Time signal
axs[0].plot(t_fft, heave_fft)
axs[0].set_title("Heave motion (last portion)")
axs[0].set_xlabel("Time [s]")
axs[0].set_ylabel("Heave [m]")
axs[0].grid()

# FFT vs frequency
axs[1].plot(freqs, amp)
axs[1].set_title("FFT amplitude spectrum")
axs[1].set_xlabel("Frequency [Hz]")
axs[1].set_ylabel("Amplitude")
axs[1].grid()

# FFT vs period
axs[2].plot(periods, amp)
axs[2].axvline(T, color='r', linestyle='--', label='Encounter period')
axs[2].set_xlim(0, 2*T)
axs[2].set_title("FFT vs period")
axs[2].set_xlabel("Period [s]")
axs[2].set_ylabel("Amplitude")
axs[2].legend()
axs[2].grid()

plt.tight_layout()


# ================= PLOTT ALL MOTION SIGNALS =================
fig, axs = plt.subplots(2, 3, figsize=(18, 8))

axs[0, 0].plot(data[:, 0], data[:, 1]/zeta)
axs[0, 0].set_title('surge')
axs[0, 1].plot(data[:, 0], data[:, 2]/zeta)
axs[0, 1].set_title('sway')
axs[0, 2].plot(data[:, 0], data[:, 3]/zeta)
axs[0, 2].set_title('heave')
axs[1, 0].plot(data[:, 0], data[:, 4]/zeta)
axs[1, 0].set_title('roll')
axs[1, 1].plot(data[:, 0], data[:, 5]/zeta)
axs[1, 1].set_title('pitch')
axs[1, 2].plot(data[:, 0], data[:, 6]/zeta)
axs[1, 2].set_title('yaw')


plt.tight_layout()
for ax in axs.flatten():
    ax.grid()

plt.show()
