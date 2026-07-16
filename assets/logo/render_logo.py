"""Render the topon logo with OVITO Tachyon (run with C:/v/ovito/Scripts/python.exe).

Follows the ovito-render skill's technical patterns (per-particle Color +
Radius via create_property; Tachyon + ambient occlusion; white background).

DELIBERATE deviation from the gallery house camera: a logo's text has to read,
and the house 3/4 view (camera_dir=(-0.55,-0.75,-0.45)) shears the letterforms
into illegibility. Renders a front-on variant and a gently tilted one.

Types (from make_logo_data.py):
    1 quiet strand   2 accent strand   3 quiet junction   4 accent junction
"""
import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

DATA = "topon_logo.data"

QUIET_STRAND = (0.788, 0.808, 0.831)     # #c9ced4
QUIET_JUNCT  = (0.682, 0.714, 0.749)     # #aeb6bf
ACCENT_STRAND = (0.055, 0.486, 0.420)    # #0e7c6b teal
ACCENT_JUNCT  = (0.043, 0.365, 0.314)    # #0b5d50

R_STRAND, R_JUNCT = 0.075, 0.17
BOND_W = 0.05


def make_pipeline():
    pipe = import_file(DATA, atom_style="full")

    def paint(frame, data):
        t = data.particles["Particle Type"].array
        colors = np.zeros((data.particles.count, 3))
        colors[t == 1] = QUIET_STRAND
        colors[t == 2] = ACCENT_STRAND
        colors[t == 3] = QUIET_JUNCT
        colors[t == 4] = ACCENT_JUNCT
        data.particles_.create_property("Color", data=colors)

        radii = np.full(data.particles.count, R_STRAND)
        radii[(t == 3) | (t == 4)] = R_JUNCT
        data.particles_.create_property("Radius", data=radii)

    pipe.modifiers.append(paint)
    pipe.add_to_scene()
    data = pipe.compute()
    if data.particles.bonds is not None:
        data.particles.bonds.vis.width = BOND_W
        data.particles.bonds.vis.use_particle_colors = True
    return pipe


def render(pipe, out, cam_dir, cam_type=Viewport.Type.Ortho, size=(1800, 620)):
    vp = Viewport(type=cam_type)
    vp.camera_dir = cam_dir
    vp.zoom_all(size=size)
    r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=16,
                        shadows=True, direct_light_intensity=1.1)
    vp.render_image(filename=out, size=size, background=(1, 1, 1), renderer=r)
    print(f"[OK] {out}")


if __name__ == "__main__":
    p = make_pipeline()
    render(p, "logo_front.png", (0.0, 0.0, -1.0))                  # text reads flat
    render(p, "logo_tilt.png", (-0.22, -0.30, -1.0))               # slight depth
    print("done")
