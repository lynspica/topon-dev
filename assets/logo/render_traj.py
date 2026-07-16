"""Render the KG trajectory to frames + assemble the animation.
Run with C:/v/ovito/Scripts/python.exe

Colours come from the SOURCE .data types (the dump only carries type ids, which
we preserve), so the wordmark stays teal throughout while the physics treats
every bead identically.
"""
import os
import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

SIZE = (1000, 340)          # per-frame; small = fast, GIF-appropriate
EVERY = 2                   # render every Nth frame (111 -> ~56)

QUIET_S = (0.788, 0.808, 0.831)
QUIET_J = (0.682, 0.714, 0.749)
ACC_S   = (0.055, 0.486, 0.420)
ACC_J   = (0.043, 0.365, 0.314)

os.makedirs("frames", exist_ok=True)
pipe = import_file("traj.lammpstrj")
n = pipe.source.num_frames
print(f"trajectory: {n} frames", flush=True)


def paint(frame, data):
    t = data.particles["Particle Type"].array
    c = np.zeros((data.particles.count, 3)); r = np.full(data.particles.count, 0.36)
    c[t == 1] = QUIET_S
    c[t == 2] = ACC_S
    c[t == 3] = QUIET_J; r[t == 3] = 0.80
    c[t == 4] = ACC_J;   r[t == 4] = 0.80
    data.particles_.create_property("Color", data=c)
    data.particles_.create_property("Radius", data=r)


pipe.modifiers.append(paint)
pipe.add_to_scene()

vp = Viewport(type=Viewport.Type.Ortho)
vp.camera_dir = (0.0, 0.0, -1.0)
vp.camera_up = (0.0, 1.0, 0.0)
vp.zoom_all(size=SIZE)          # lock the camera on frame 0 -> no jitter
r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=8, shadows=False)

for i in range(0, n, EVERY):
    out = f"frames/f{i:04d}.png"
    vp.render_image(filename=out, size=SIZE, background=(1, 1, 1),
                    renderer=r, frame=i)
    print(f"[{i}/{n}] {out}", flush=True)
print("frames done")
