import numpy as np
import os

# =============================================================================
# 1. PARAMETERS & INITIALIZATION
# =============================================================================
heading  = 0.0
head_ang = heading * np.pi / 180.0

rho   = 1000.0
zeta  = 0.025
g     = 9.81
lam   = 3.0
k     = 2 * np.pi / lam
Ucur  = 0.2

omega = np.sqrt(g * k) + Ucur * k * np.cos(head_ang)
T_e   = 2 * np.pi / omega  # Encounter period

Lpp   = 1.0
B     = 1.0

# Non-dimensionalization denominators
denom    = rho * g * zeta**2 * B**2 *0.5
denomMom = denom * Lpp

print(f"lam/Lpp = {lam/Lpp:.4f}")
print(f"tau = {omega*Ucur/9.81:.4f}")

# =============================================================================
# 2. LOAD OPENFOAM DATA
# =============================================================================
fnameForce  = "postProcessing/MWLforcesClaude/0/force.dat"
fnameMoment = "postProcessing/MWLforcesClaude/0/moment.dat"

# Handle potential missing files gracefully
try:
    data    = np.loadtxt(fnameForce, skiprows=1)
    dataMom = np.loadtxt(fnameMoment, skiprows=1)
except FileNotFoundError as e:
    print(f"Error: Could not find force/moment data files. {e}")
    exit()

time = data[:, 0]
Fx   = data[:, 1]
Fy   = data[:, 2]
Mz   = dataMom[:, 3]

# =============================================================================
# 3. HARMONIC LEAST SQUARES FIT (Replaces Zero-Crossings)
# =============================================================================
def extract_mean_harmonic(t, signal, omega_e, num_harmonics=2):
    """
    Fits a Fourier series to the signal to exactly extract the mean (0th harmonic).
    Equation: F(t) = A0 + A1*cos(w*t) + B1*sin(w*t) + A2*cos(2*w*t) + ...
    """
    M = np.ones((len(t), 1))
    for n in range(1, num_harmonics + 1):
        M = np.column_stack((M, np.cos(n * omega_e * t), np.sin(n * omega_e * t)))
        
    coeffs, _, _, _ = np.linalg.lstsq(M, signal, rcond=None)
    return coeffs[0] # Return A0 (the exact mean drift load)

# Define time window for steady-state evaluation (e.g., last 4 periods)
num_steady_periods = 4
eval_time_window = num_steady_periods * T_e

# Slice the data to only include the steady-state tail
t_end = time[-1]
t_start = max(time[0], t_end - eval_time_window) # Fallback if simulation is short
steady_idx = np.searchsorted(time, t_start)

t_steady  = time[steady_idx:]
Fx_steady = Fx[steady_idx:]
Fy_steady = Fy[steady_idx:]
Mz_steady = Mz[steady_idx:]

# Extract raw global means
mean_Fx = extract_mean_harmonic(t_steady, Fx_steady, omega)
mean_Fy = extract_mean_harmonic(t_steady, Fy_steady, omega)
mean_Mz = extract_mean_harmonic(t_steady, Mz_steady, omega)

# =============================================================================
# 4. ROTATION & NON-DIMENSIONALIZATION
# =============================================================================
# Project forces into the ship's local coordinate system based on heading
meanFx_p = mean_Fx * np.cos(head_ang) - mean_Fy * np.sin(head_ang)
meanFy_p = mean_Fx * np.sin(head_ang) + mean_Fy * np.cos(head_ang)
meanMz_p = mean_Mz # Z-axis moment remains unchanged under XY rotation

# Non-dimensionalize
meanFx_p_ND = meanFx_p / denom
meanFy_p_ND = meanFy_p / denom
meanMz_ND   = meanMz_p / denomMom

# =============================================================================
# 5. NICE PRINT-OUTS
# =============================================================================
print("\n" + "="*60)
print(" 2ND-ORDER MEAN WAVE LOADS (FOURIER FIT)")
print("="*60)
print(f" Evaluated over last {num_steady_periods} periods")
print(f" Time Window : {t_start:.2f} s to {t_end:.2f} s")
print(f" Enc. Period : {T_e:.4f} s")
print("-" * 60)
print(" [1] RAW MEAN FORCES (Global Frame)")
print(f"     Fx : {mean_Fx:>10.4f} N")
print(f"     Fy : {mean_Fy:>10.4f} N")
print(f"     Mz : {mean_Mz:>10.4f} Nm")
print("-" * 60)
print(f" [2] ROTATED MEAN FORCES (Heading = {heading} deg)")
print(f"     Fx': {meanFx_p:>10.4f} N")
print(f"     Fy': {meanFy_p:>10.4f} N")
print(f"     Mz': {meanMz_p:>10.4f} Nm")
print("-" * 60)
print(" [3] NON-DIMENSIONAL LOADS")
print(f"     Fx'/ND : {meanFx_p_ND:>10.5f}")
print(f"     Fy'/ND : {meanFy_p_ND:>10.5f}")
print(f"     Mz'/ND : {meanMz_ND:>10.5f}")
print("="*60 + "\n")

# =============================================================================
# 6. FILE APPENDING & PLOTTING
# =============================================================================
with open('../res30U.dat', 'a') as f:
    f.write(f"{lam/Lpp} {meanFx_p_ND} {meanFy_p_ND} {meanMz_ND}\n")

# Load experimental/Wadam data
expFolder = os.path.expanduser('~/expWadRes/')
try:
    wamDat = np.loadtxt(os.path.join(expFolder, 'wadDat30.dat'))
    expFx  = np.loadtxt(os.path.join(expFolder, 'expFx30U0.dat'))
    expFy  = np.loadtxt(os.path.join(expFolder, 'expFy30U0.dat'))
    expMz  = np.loadtxt(os.path.join(expFolder, 'expMz30U0.dat'))
    expFxU = np.loadtxt(os.path.join(expFolder, 'expFx30U0387.dat'))
    expFyU = np.loadtxt(os.path.join(expFolder, 'expFy30U0387.dat'))
    expMzU = np.loadtxt(os.path.join(expFolder, 'expMz30U0387.dat'))
except FileNotFoundError as e:
    print(f"Warning: Could not load experimental comparison data. {e}")
    # Create empty arrays to prevent plotting crashes if files are missing
    wamDat = expFx = expFy = expMz = expFxU = expFyU = expMzU = np.zeros((1,4))


