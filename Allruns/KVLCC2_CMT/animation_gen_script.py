"""
ParaView orbit-camera animation driven by physical simulation time.

Camera azimuth:
    theta(t) = theta0                                  for t <  t_start
    theta(t) = theta0 + omega * (t - t_start)          for t >= t_start

Run in the ParaView Python Shell (with your pipeline already loaded),
or with pvpython after loading your state file (see LOAD_STATE below).

Encode frames afterwards, e.g.:
    ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p orbit.mp4
"""

import os
import numpy as np
from paraview.simple import *

# ----------------------------------------------------------------------
# User parameters
# ----------------------------------------------------------------------
LOAD_STATE = None            # e.g. 'case.pvsm' when running with pvpython; None in GUI shell

t_start = 8.37               # [s] simulation time at which rotation begins
omega   = 0.1               # [rad/s] azimuthal rate, in simulation seconds (sign = direction)
theta0  = -np.pi / 2        # [rad] initial azimuth (measured from +x about the axis)

focal   = np.array([11.2, 0.0, 0.0])   # [m] point the camera looks at
axis    = np.array([0.0, 0.0, 1.0])    # orbit axis (also used as view-up)
R       = 40.0                        # [m] orbit radius
h       = 20.0                         # [m] camera height along axis above focal point

# Frame sampling:
#   None  -> one frame per stored data time step
#   float -> uniform frame spacing in simulation time [s]; data snaps to the
#            latest available time step (use for non-uniform write intervals)
dt_frame = None
t_end    = None              # [s] stop time; None -> last data time step

out_dir    = r'C:\Users\nikospe\Documents\orbit_frames'
resolution = [1920, 1080]

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

t_first = data_times[0]
t_last  = data_times[-1] if t_end is None else min(t_end, data_times[-1])

if dt_frame is None:
    frame_times = data_times[data_times <= t_last + 1e-12]
else:
    frame_times = np.arange(t_first, t_last + 0.5 * dt_frame, dt_frame)

# Orthonormal basis (e1, e2) in the plane normal to the axis
a = axis / np.linalg.norm(axis)
ref = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
e1 = ref - np.dot(ref, a) * a
e1 /= np.linalg.norm(e1)
e2 = np.cross(a, e1)

def azimuth(t):
    return theta0 + omega * max(0.0, t - t_start)

def data_time_at(t):
    """Latest stored data time <= t (ParaView would otherwise snap similarly)."""
    idx = np.searchsorted(data_times, t + 1e-12) - 1
    return data_times[max(idx, 0)]

os.makedirs(out_dir, exist_ok=True)

# ----------------------------------------------------------------------
# Render loop
# ----------------------------------------------------------------------
print(f'{len(frame_times)} frames, t = {frame_times[0]:.3f} .. {frame_times[-1]:.3f} s')
print(f'Rotation from t = {t_start} s at {omega} rad/s '
      f'({np.degrees(omega):.2f} deg/s, period {2*np.pi/abs(omega) if omega else np.inf:.1f} s)')

for i, t in enumerate(frame_times):
    scene.AnimationTime = float(data_time_at(t))

    th  = azimuth(t)
    pos = focal + R * (np.cos(th) * e1 + np.sin(th) * e2) + h * a

    view.CameraPosition   = pos.tolist()
    view.CameraFocalPoint = focal.tolist()
    view.CameraViewUp     = a.tolist()
    view.GetRenderer().ResetCameraClippingRange()

    Render(view)
    SaveScreenshot(os.path.join(out_dir, f'frame_{i:04d}.png'), view,
                   ImageResolution=resolution)

    if i % 20 == 0:
        print(f'  frame {i:4d}  t = {t:8.3f} s  theta = {np.degrees(th):8.2f} deg')

print('Done.')