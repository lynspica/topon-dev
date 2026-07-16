"""Render the topon logo variants + a labelled 2x2 comparison sheet.
Run with C:/v/ovito/Scripts/python.exe

Follows the ovito-render skill's technical patterns (per-particle Color/Radius
via create_property, Tachyon + ambient occlusion, white background). Deliberate
deviation from the gallery house camera: a wordmark has to READ, and the house
3/4 view shears the letterforms illegible -- so this is front-on.
"""
import sys
import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

VARIANTS = ["clean", "copolymer", "entangled", "organic"]
SIZE = (1300, 440)          # comparison size; final pick re-rendered bigger

COL = {
    1: (0.788, 0.808, 0.831),   # quiet strand   #c9ced4
    3: (0.682, 0.714, 0.749),   # quiet junction #aeb6bf
    2: (0.055, 0.486, 0.420),   # accent-A strand   teal #0e7c6b
    4: (0.043, 0.365, 0.314),   # accent-A junction teal #0b5d50
    5: (0.878, 0.478, 0.325),   # accent-B strand   coral #e07a53
    6: (0.776, 0.373, 0.227),   # accent-B junction coral #c65f3a
    7: (0.851, 0.467, 0.024),   # entanglement ring amber #d97706
}
RAD = {1: 0.075, 2: 0.075, 5: 0.075, 3: 0.17, 4: 0.17, 6: 0.17, 7: 0.13}


def render(variant):
    pipe = import_file(f"topon_{variant}.data", atom_style="full")

    def paint(frame, data):
        t = data.particles["Particle Type"].array
        colors = np.zeros((data.particles.count, 3))
        radii = np.full(data.particles.count, 0.075)
        for tid, c in COL.items():
            m = t == tid
            if m.any():
                colors[m] = c
                radii[m] = RAD[tid]
        data.particles_.create_property("Color", data=colors)
        data.particles_.create_property("Radius", data=radii)

    pipe.modifiers.append(paint)
    pipe.add_to_scene()
    d = pipe.compute()
    if d.particles.bonds is not None:
        d.particles.bonds.vis.width = 0.05
        d.particles.bonds.vis.use_particle_colors = True

    vp = Viewport(type=Viewport.Type.Ortho)
    vp.camera_dir = (0.0, 0.0, -1.0)
    vp.zoom_all(size=SIZE)
    r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_samples=14,
                        shadows=True, direct_light_intensity=1.1)
    out = f"var_{variant}.png"
    vp.render_image(filename=out, size=SIZE, background=(1, 1, 1), renderer=r)
    pipe.remove_from_scene()
    print(f"[OK] {out}", flush=True)
    return out


if __name__ == "__main__":
    want = sys.argv[1:] or VARIANTS
    for v in want:
        render(v)
    print("renders done")
