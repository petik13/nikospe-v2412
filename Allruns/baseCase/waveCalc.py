import numpy as np
import matplotlib.pyplot as plt

lam = 2.0
steepness = 1/60
Ucur = 0.0
heading = 0.0

head_ang = np.radians(heading)
k = lam / (2.0 * np.pi)
g = 9.81
omega = np.sqrt(g * k) + Ucur * k * np.cos(head_ang)
T = 2.0 * np.pi / omega

zeta = 0.5*steepness*lam


print(f'Wave period: {T:.2f} s')
print(f'Wave amplitude: {zeta:.2f} m')
