"""
ParaView: render body-fixed CFD results as seen by an EARTH-FIXED camera.

Kinematics (planar, yaw only, z up, positive r = counter-clockwise seen from above):
  Reference point C (body coords x_c) moves with body-frame velocity (u, v),
  body yaws at rate r, starting at t_start. Before t_start: straight course.

  tau = t - t_start
  psi(t) = r * max(tau, 0)
  p_C(t) = x_c + int_0^tau R(psi(s)) [u, v] ds
         = x_c + [u*tau, v*tau]                                         tau <= 0
         = x_c + (1/r) [ u sin(r tau) - v (1 - cos(r tau)),
                         u (1 - cos(r tau)) + v sin(r tau) ]           tau  > 0

  Earth frame is chosen to coincide with the body frame at t = t_start, so
  camera / focal point are specified in body coordinates as they are at t_start.

  Earth -> body transform (what the camera must be set to in the CFD frame):
  x_b = x_c + R(psi)^T (x_E - p_C(t))

Encode afterwards:
  ffmpeg -framerate 30 -i frame_%04d.png -c:v libx264 -crf 18 -pix_fmt yuv420p earth_view.mp4
"""

import os
import numpy as np
from paraview.simple import *

# ----------------------------------------------------------------------
# User parameters
# ----------------------------------------------------------------------
LOAD_STATE = None

# Body motion
t_start = 8.37                      # [s] yaw rotation starts here
u       = -0.33                      # [m/s] surge velocity of point C (body x)
v       = 0.0                       # [m/s] sway velocity of point C (body y), e.g. drift
r       = -0.1                       # [rad/s] yaw rate (positive = turn to +y / port, z up)
x_c     = np.array([11.2, 0.0, 0.0]) # [m] body-frame point whose velocity is (u, v)
                                    #     = rotation centre of the CFD frame (SRF origin / CoG)

# Earth-fixed camera, given in body coordinates AT t_start
if abs(r) > 1e-12:
    turn_centre = x_c + np.array([-v / r, u / r, 0.0])   # centre of the steady turning circle
else:
    turn_centre = x_c.copy()

focal_E = turn_centre.copy()                              # [m] what the camera looks at
cam_E   = focal_E + np.array([0.0, -14.0, 20.0])          # [m] camera position
up_E    = np.array([0.0, 0.0, 1.0])                       # view-up (must not be parallel to view dir)

# Frame sampling
dt_frame = None            # None -> one frame per data step; float -> uniform spacing [s]
t_end    = None

out_dir    = r'C:\Users\nikospe\Documents\earth_frames'
resolution = [1920, 1080]

# ----------------------------------------------------------------------
# Kinematics
# ----------------------------------------------------------------------
def body_pose(t):
    """Heading psi [rad] and earth position of C (3-vector) at time t."""
    tau = t - t_start
    if tau <= 0.0:
        return 0.0, x_c + np.array([u * tau, v * tau, 0.0])
    if abs(r) < 1e-12:
        return 0.0, x_c + np.array([u * tau, v * tau, 0.0])
    s, c = np.sin(r * tau), np.cos(r * tau)
    dp = np.array([(u * s - v * (1.0 - c)) / r,
                   (u * (1.0 - c) + v * s) / r,
                   0.0])
    return r * tau, x_c + dp

def earth_to_body(x_E, psi, p_C):
    """Point in earth frame -> body (CFD) frame."""
    d = x_E - p_C
    c, s = np.cos(psi), np.sin(psi)
    return x_c + np.array([ c * d[0] + s * d[1],
                           -s * d[0] + c * d[1],
                            d[2]])

def earth_dir_to_body(d_E, psi):
    """Direction vector earth -> body (rotation only)."""
    c, s = np.cos(psi), np.sin(psi)
    return np.array([c * d_E[0] + s * d_E[1], -s * d_E[0] + c * d_E[1], d_E[2]])

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------
if LOAD_STATE:
    LoadState(LOAD_STATE)

view  = GetActiveViewOrCreate('RenderView')
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
data_times = np.array(GetTimeKeeper().TimestepValues, dtype=float)
if data_times.size == 0:
    raise RuntimeError('No time steps found in the pipeline.')

t_last = data_times[-1] if t_end is None else min(t_end, data_times[-1])
if dt_frame is None:
    frame_times = data_times[data_times <= t_last + 1e-12]
else:
    frame_times = np.arange(data_times[0], t_last + 0.5 * dt_frame, dt_frame)

def data_time_at(t):
    idx = np.searchsorted(data_times, t + 1e-12) - 1
    return data_times[max(idx, 0)]

os.makedirs(out_dir, exist_ok=True)

if abs(r) > 1e-12:
    print(f'Turning radius of C: {np.hypot(u, v) / abs(r):.3f} m, '
          f'period {2 * np.pi / abs(r):.1f} s, {np.degrees(r):.2f} deg/s')
print(f'{len(frame_times)} frames, t = {frame_times[0]:.3f} .. {frame_times[-1]:.3f} s')

# ----------------------------------------------------------------------
# Render loop
# ----------------------------------------------------------------------
for i, t in enumerate(frame_times):
    scene.AnimationTime = float(data_time_at(t))

    psi, p_C = body_pose(t)
    view.CameraPosition   = earth_to_body(cam_E,   psi, p_C).tolist()
    view.CameraFocalPoint = earth_to_body(focal_E, psi, p_C).tolist()
    view.CameraViewUp     = earth_dir_to_body(up_E, psi).tolist()

    Render(view)
    SaveScreenshot(os.path.join(out_dir, f'frame_{i:04d}.png'), view,
                   ImageResolution=resolution)

    if i % 20 == 0:
        print(f'  frame {i:4d}  t = {t:8.3f} s  psi = {np.degrees(psi):8.2f} deg  '
              f'C_E = ({p_C[0]:7.3f}, {p_C[1]:7.3f})')

print('Done.')