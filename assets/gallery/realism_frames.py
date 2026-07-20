"""'Making the junctions realistic' — SC -> +BCC/FCC -> +jitter.

Uses the network-lattice-realism skill's own point-set and Gaussian-matching code
(imported from the skill's scripts, not reimplemented), so the numbers on screen
are the numbers that skill reports.

The story, in the skill's words: a perfect lattice has CRYSTALLINE junctions;
mixing SC/BCC/FCC adds neighbour shells but *by itself leaves them crystalline*;
site jitter is the primary realism lever that makes them liquid-like. The
animation shows exactly that, with a live junction g(r).

Stages:
  A  SC only                      g(r) peak ~6.8  crystalline
  B  + BCC + FCC sublattices      g(r) peak ~5.0  still crystalline
  C  + site jitter ramped to 0.25 g(r) peak ~1.05 liquid-like

TOPON-python stage: writes f{k}.data (junctions as atoms, Gaussian-matched edges
as bonds) + meta.json (per-frame g(r) curve, peak, S(q->0), stage label).
"""
import json
import sys
from pathlib import Path

import numpy as np

SKILL = Path(r"C:/Users/ahmet/OneDrive - Northwestern University/DOW-Ahmet/Skills"
             r"/network-lattice-realism/scripts")
sys.path.insert(0, str(SKILL))
from build_and_check import (lattice, make_pointset, add_jitter, bmatch,
                             junction_gr_peak, Sq0)

M = 6                 # M^3 cells
MIX = (0.2, 0.4, 0.4)  # SC / BCC / FCC fractions (skill's suggestion at N=80)
JITTER = 0.25
RMS_LAT = 1.75        # target rms Ree in spacing units (skill: N=80, rho=0.30)
CUTOFF, EXT_CAP, QUOTA = 2.6, 3.2, 4


def gr_curve(pos, L, rmax=2.2, nbins=44):
    """Full radial distribution function of the junction point set."""
    n = len(pos); rho = n / L**3
    edges = np.linspace(0.25, rmax, nbins + 1)
    d2 = []
    for i in range(n):
        d = pos[i+1:] - pos[i]; d -= L * np.round(d / L)
        d2.append((d**2).sum(1))
    d = np.sqrt(np.concatenate(d2))
    h, _ = np.histogram(d, bins=edges)
    shell = 4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
    g = h / (0.5 * n * rho * shell)
    return 0.5 * (edges[:-1] + edges[1:]), g


def write_data(path, pos, L, edges):
    N = len(pos)
    lines = ["realism frame\n", f"{N} atoms", f"{len(edges)} bonds",
             "1 atom types", "1 bond types",
             f"0 {L:.4f} xlo xhi", f"0 {L:.4f} ylo yhi", f"0 {L:.4f} zlo zhi",
             "", "Masses", "", "1 1.0", "", "Atoms # full", ""]
    for i, (x, y, z) in enumerate(pos):
        lines.append(f"{i+1} 1 1 0 {x:.5f} {y:.5f} {z:.5f}")
    lines += ["", "Bonds", ""]
    for k, (a, b) in enumerate(edges):
        lines.append(f"{k+1} 1 {int(a)+1} {int(b)+1}")
    Path(path).write_text("\n".join(lines) + "\n")


def build_all(outdir, n_jitter_frames=20, hold=6, seed=12345):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    # --- the two point sets (same box, so the camera never moves) ---
    P_sc, L = make_pointset("sc", M, MIX, np.random.default_rng(seed))
    P_mix, L_mix = make_pointset("mix", M, MIX, np.random.default_rng(seed))
    # normalise the mix box to the SC box so both render at one scale
    P_mix = P_mix * (L / L_mix)

    frames = []          # (label, points, jitter_sigma)
    for _ in range(hold):
        frames.append(("A", P_sc, 0.0))
    for _ in range(hold):
        frames.append(("B", P_mix, 0.0))
    # Ramp starts ABOVE zero: linspace(0, ...) would make the first stage-C
    # frame identical to stage B (jitter 0) but labelled "add site jitter",
    # which reads as a mislabelled duplicate.
    for j in np.linspace(0, JITTER, n_jitter_frames + 1)[1:]:
        frames.append(("C", P_mix, float(j)))
    for _ in range(hold):
        frames.append(("C", P_mix, JITTER))

    meta = []
    for k, (stage, P, sig) in enumerate(frames):
        # one jitter draw per sigma, from a fixed seed, so the ramp is smooth
        # (each site walks steadily outward instead of re-randomising each frame)
        r2 = np.random.default_rng(999)
        step = r2.normal(size=P.shape)
        pos = (P + sig * step) % L if sig > 0 else P.copy()
        e, _d = bmatch(pos, L, QUOTA, RMS_LAT, CUTOFF, EXT_CAP,
                       np.random.default_rng(7))
        r, g = gr_curve(pos, L)
        peak = junction_gr_peak(pos, L)
        sq = Sq0(pos, L)
        write_data(outdir / f"f{k:03d}.data", pos, L, e)
        meta.append({"stage": stage, "jitter": sig, "n": int(len(pos)),
                     "edges": int(len(e)), "peak": float(peak), "sq0": float(sq),
                     "r": r.tolist(), "g": [float(x) for x in g]})
        print(f"  f{k:03d} {stage} jitter {sig:.3f}: {len(pos)} junctions, "
              f"{len(e)} edges, g(r) peak {peak:.2f}, S(0) {sq:.3f}", flush=True)

    (outdir / "meta.json").write_text(json.dumps(
        {"L": float(L), "M": M, "mix": list(MIX), "jitter": JITTER,
         "frames": meta}, indent=1))
    print(f"[realism_frames] {len(frames)} frames -> {outdir}")


if __name__ == "__main__":
    build_all(sys.argv[1] if len(sys.argv) > 1 else "C:/tmp/realism_frames")
