"""Generate the systems the gallery renders, straight from topon's Python API.

Every panel sits on a STRICT-SCULPTED, heterogeneous network -- junction
functionality 2-6, mean 4.0 -- not a perfect lattice. topon's sculpting stage
thins a full SC lattice to a target edge count (`degree_distribution = "e:N"`),
which is what gives the real degree spread.

Four systems:
  sculpt_250         5x5x5 SC, 375 -> 250 edges          -- the CG hero arc
  atom_sculpt        3x3x3 PDMS, 81 -> 54 edges          -- the atomistic arc
  copoly_{block,      4x4x4, 192 -> 128 edges, A/B 50:50 -- sequence panels
    random,alternating}
  entangled_grafted  5x5x5, -> 250 edges, 12 entanglements + grafts

The two arcs (sculpt_250, atom_sculpt) then go through LAMMPS; the copolymer and
entanglement panels use the 03_Conformation (lattice) state only.

Systems are big and transient, so they are written OUTSIDE the repo. Point
TOPON_GALLERY_SYSTEMS at a scratch dir (defaults to ./systems next to this file).

Usage:  python gen_systems.py [name ...]      # default: all four
"""
import json
import os
import sys
from pathlib import Path

from topon.config import load_config
from topon.pipeline import Pipeline

OUT = Path(os.environ.get("TOPON_GALLERY_SYSTEMS", Path(__file__).parent / "systems"))
OUT.mkdir(parents=True, exist_ok=True)


def cg(lattice, edges, assignment=None):
    d = {
        "study": {"name": None, "output_dir": None},
        "topology": {"source": "generate", "generator": {
            "lattice_size": lattice, "lattice_type": "SC", "max_functionality": 6,
            "degree_distribution": f"e:{edges}"}},
        "chemistry": {"model_type": "coarse_grained", "bead_density": 0.85,
                      "degree_of_polymerization": 20},
        "assignment": assignment or {},
        "conformation": {"overlap_cutoff": 0.01, "overlap_max_iters": 10},
        "simulation": {"include_angles": True, "pair_style": "attractive"},
        "execution": {"auto_run": False},
    }
    return d


def copoly(arrangement):
    return cg("4x4x4", 128, {"copolymer": {"enabled": True, "per_edge_type": {"A": {
        "arrangement": arrangement,
        "composition": [{"monomer": "A", "fraction": 0.5},
                        {"monomer": "B", "fraction": 0.5}]}}}})


CONFIGS = {
    # CG hero: 5x5x5 sculpted 375 -> 250 (mean functionality 4.0).
    "sculpt_250": cg("5x5x5", 250),
    # Copolymer sequences on a sculpted 4x4x4 (192 -> 128, mean f 4.0).
    "copoly_block":       copoly("block"),
    "copoly_random":      copoly("random"),
    "copoly_alternating": copoly("alternating"),
    # Entanglements + grafts on the same 5x5x5 sculpt as the hero. NOTE: the
    # closest crossing varies run to run (partner selection is unseeded); the
    # shipped realisation was picked from six for a tight 0.39 sigma pair.
    "entangled_grafted": cg("5x5x5", 250, {
        "entanglements": {"enabled": True, "target": 12, "target_type": "count"},
        "grafts": {"enabled": True, "per_edge_type": {
            "A": {"graft_density": 0.05, "side_chain_dp": 5,
                  "side_chain_monomer": "G"}}}}),
    # Atomistic arc: 3x3x3 PDMS sculpted 81 -> 54. NOTE e:54, not e:128 -- a
    # 3x3x3 SC has only 81 edges, so a larger target never terminates.
    "atom_sculpt": {
        "study": {"name": None, "output_dir": None},
        "topology": {"source": "generate", "generator": {
            "lattice_size": "3x3x3", "lattice_type": "SC", "max_functionality": 6,
            "degree_distribution": "e:54"}},
        "chemistry": {"model_type": "atomistic", "monomer_type": "PDMS",
                      "degree_of_polymerization": 20, "density": 0.85},
        "conformation": {"overlap_cutoff": 0.01, "overlap_max_iters": 10},
        "simulation": {"include_angles": True},
        "execution": {"auto_run": False},
    },
}


def build(name):
    d = json.loads(json.dumps(CONFIGS[name]))
    d["study"]["name"] = name
    d["study"]["output_dir"] = str(OUT / name)
    cfg_path = OUT / f"{name}.json"
    cfg_path.write_text(json.dumps(d, indent=1))
    print(f"\n=== {name} ===", flush=True)
    Pipeline(load_config(str(cfg_path))).run()


if __name__ == "__main__":
    for nm in (sys.argv[1:] or list(CONFIGS)):
        try:
            build(nm)
        except Exception as e:
            print(f"[FAIL] {nm}: {type(e).__name__}: {e}", flush=True)
    print("\ngen done")
