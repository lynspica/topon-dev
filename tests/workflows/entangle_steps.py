"""Entanglement, built up one step at a time on one fixed CG network.

Every step uses the *same* lattice, the same chains, the same box and the
same LAMMPS protocol. Only the set of prescribed contacts changes. That is
the point: when a step misbehaves, the lattice is not a suspect.

    python tests/workflows/entangle_steps.py --step 1      # lattice only, no braids
    python tests/workflows/entangle_steps.py --step 1 --run-md

Steps (later ones are added as each is agreed):

    1  the network itself: mixed lattice, CG chains, linear paths, MD runs
    2  one pair of chains entangled once
    3  several entanglements along one chain
    4  one chain entangled with several different chains
    5  a small composite: A-B, A-C, D-B twice, D-E once

Output per step lands in tests/output/entangle_steps/step<N>/, laid out the
way the rest of the pipeline lays things out (02_Chemistry, 03_Conformation,
04_Simulation), so the LAMMPS scripts are the generated ones and not a
special case written for this script.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import Chem  # noqa: E402

from topon.conformation import ConformationManager  # noqa: E402
from topon.conformation.junction_shell import (  # noqa: E402
    apply_junction_shells,
)
from topon.conformation.entanglement import (  # noqa: E402
    BraidShape,
    ContactRequest,
    allocate_contacts,
    closest_approach,
    compose_chain_path,
    far_closed_linking,
    gap_at,
    min_separation,
)
from topon.simulation import SimulationRunner  # noqa: E402
from topon.topology.generator_python import PythonTopologyGenerator  # noqa: E402
from topon.utils import write_lammps_displacement_file  # noqa: E402
from topon.writers import CGWriter, LammpsInputGenerator  # noqa: E402

OUT = ROOT / "tests/output/entangle_steps"

# One network for the whole series. Chosen, not tuned: 4x4x4 is large enough
# that a chain has neighbours in every shell but small enough to look at, and
# the SC/BCC/FCC mix gives strand pairs at several different separations,
# which is exactly what "1st, 2nd, 3rd neighbour" needs to mean something.
LATTICE = dict(lattice="MIX", dims=(4, 4, 4), mix={"SC": 0.2, "BCC": 0.4, "FCC": 0.4},
               cutoff=1.0, max_func=12, degree_dist="0:0,1:0", seed=42)
DP = 40                 # beads per chain between its two junctions
DENSITY = 0.85          # melt density, for the --density comparison run

# Bond length every chain is built at, and the number that fixes the scale:
# the longest chain, drawn straight with DP beads, gets bonds this long.
#
# The scale follows from the geometry, not from a chosen density; the density
# is reported instead. 0.90 sits just under the Kremer-Grest equilibrium of
# ~0.97, so the protocol has nothing to stretch.
BOND = 0.90

# Minimum separation asked of the first beads of chains sharing a junction.
# Set to 0 to build without the shell, which is how the overlap counts in
# junction_shell's docstring were measured.
SHELL_SPACING = 1.0

# The scale is pinned by the LONGEST path in the system, not by the longest
# chord and not per chain. Only the longest path can reach the FENE limit, so
# holding it at BOND puts every other chain below BOND automatically: shorter
# path, same bead count, shorter bonds. One DP for the whole network, and no
# bond anywhere can break.
#
# Once entanglements exist the longest path is a braided chain, not a chord,
# because a braid adds contour. Measured on BraidShape(): one winding covers
# 7.19 sigma of chord with 12.50 sigma of path. Sizing from chords instead
# would leave exactly those chains over the limit -- with DP fixed at 40 and
# the scale set from the longest chord, braids came out needing 4.05 sigma
# bonds against a limit of 1.5.


class _Cfg:
    def __init__(self, spec):
        self.lattice_type = spec["lattice"]
        self.lattice_size = spec["dims"]
        self.max_functionality = spec["max_func"]
        self.degree_distribution = spec["degree_dist"]
        self.periodicity = "111"
        self.mix_fractions = spec["mix"]
        self.mix_cutoff = spec["cutoff"]


def build_network(spec=LATTICE):
    """The shared network. Deterministic, so every step gets the same one."""
    random.seed(spec["seed"])
    np.random.seed(spec["seed"])
    gen = PythonTopologyGenerator(_Cfg(spec))
    graphs = gen.generate(trials=4000, max_saves=1, time_limit=120)
    if not graphs:
        raise SystemExit(f"sculpting produced no {spec['lattice']} network")
    return graphs[0]


def geometry(graph, dp=DP, density=None, bond=BOND, scale=None):
    """Junction positions and chain chords in sigma, plus the lattice scale.

    Two ways to fix the scale, the same equation solved for different
    unknowns:

    ``density``  the melt route. Volume follows from bead count, and the
                 bond length is whatever falls out -- 0.150 sigma here,
                 which is a chain collapsed onto its own chord.

    ``bond``     the geometric route. The scale is set so the longest chain,
                 drawn straight with ``dp`` beads, has bonds this long.
                 Density is reported rather than chosen.

    Chords use the minimum image, so a chain whose junctions sit on
    opposite faces is the short one that crosses the boundary, not a line
    straight across the system.
    """
    box = np.asarray(graph.graph["box"], float)
    n_beads = graph.number_of_edges() * dp + graph.number_of_nodes()
    cells = float(np.prod(box))

    raw = {n: np.asarray(d["pos"], float) for n, d in graph.nodes(data=True)}
    c_max = max(
        np.linalg.norm((raw[v] - raw[u]) - box * np.round((raw[v] - raw[u]) / box))
        for u, v in graph.edges())

    if scale is None:
        if density is not None:
            scale = ((n_beads / density) / cells) ** (1.0 / 3.0)
        else:
            scale = bond * (dp + 1) / c_max
    density_out = n_beads / (cells * scale ** 3)

    pos = {n: p * scale for n, p in raw.items()}
    L = box * scale

    chords, ends = {}, {}
    for k, (u, v) in enumerate(sorted(graph.edges())):
        a = pos[u]
        mic = (pos[v] - a) - L * np.round((pos[v] - a) / L)
        chords[k] = (a, a + mic)
        ends[k] = (u, v)
    return dict(box=box, scale=scale, L=L, pos=pos, chords=chords, ends=ends,
                density=density_out, c_max=c_max * scale,
                bond=c_max * scale / (dp + 1))


def path_length(p):
    return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())


def resample(p, n):
    """Re-space a polyline onto ``n`` points at equal arc length.

    Bead count follows the path, so a detour has to be walked at the same
    bond length as a straight run rather than covered by stretching.
    """
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] < 1e-12:
        return np.repeat(p[:1], n, axis=0)
    want = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(want, s, p[:, d]) for d in range(3)])


def scale_for_longest(paths_unit, dp=DP, bond=BOND):
    """Scale in sigma per lattice unit that puts the longest path at ``bond``.

    ``paths_unit`` are the paths in lattice units. The longest of them is the
    only one that can reach the FENE limit, so it sets the scale and every
    other chain lands below ``bond`` on its own.
    """
    longest = max(path_length(p) for p in paths_unit.values())
    return bond * (dp + 1) / longest


def place_beads(paths, dp=DP, ends=None, shell=SHELL_SPACING):
    """Re-space every path onto the same dp+2 beads, ends included.

    With ``ends`` given, the beads nearest each junction are also seated on
    a spread shell, so the chains meeting there do not start on top of one
    another. See topon.conformation.junction_shell for why that decides
    whether a prescribed entanglement survives the protocol.
    """
    out = {k: resample(p, dp + 2) for k, p in paths.items()}
    if ends is not None and shell:
        out = apply_junction_shells(out, ends, spacing=shell)
    return out


def linear_paths_unit(geo):
    """Straight chain paths in lattice units: just the two endpoints."""
    s = geo["scale"]
    return {k: np.stack([a0 / s, a1 / s]) for k, (a0, a1) in geo["chords"].items()}


# ---------------------------------------------------------------------------
# Chemistry and coordinates
# ---------------------------------------------------------------------------

def write_system(graph, geo, paths, root):
    """Build the CG molecule and hand the bead coordinates to the pipeline.

    Each chain gets exactly as many interior beads as its path has interior
    points, so the chemistry follows the geometry rather than the reverse.

    Coordinates go out in lattice units and are scaled back up by the
    conformation stage, which is how the rest of the pipeline does it; that
    stage also owns the box and the wrapping, so nothing here has to.
    """
    chem = root / "02_Chemistry"
    chem.mkdir(parents=True, exist_ok=True)

    mol = Chem.RWMol()
    node_atom = {}
    for node in sorted(graph.nodes()):
        idx = mol.AddAtom(Chem.Atom("Si"))
        mol.GetAtomWithIdx(idx).SetProp("bead_type", "J")
        node_atom[node] = idx

    chain_atoms = {}
    for k, (u, v) in sorted(geo["ends"].items()):
        prev, atoms = node_atom[u], []
        for _ in range(len(paths[k]) - 2):
            idx = mol.AddAtom(Chem.Atom("Si"))
            mol.GetAtomWithIdx(idx).SetProp("bead_type", "A")
            mol.AddBond(prev, idx, Chem.BondType.SINGLE)
            atoms.append(idx)
            prev = idx
        mol.AddBond(prev, node_atom[v], Chem.BondType.SINGLE)
        chain_atoms[k] = atoms

    CGWriter(mol, str(chem / "system.data"), include_angles=True).write()
    (chem / "system.in.settings").write_text("# CG settings\n")

    s = geo["scale"]
    node_xyz = {node_atom[n]: tuple(p / s) for n, p in geo["pos"].items()}
    write_lammps_displacement_file(node_xyz, s, s, s,
                                   str(chem / "system_nodes.displace"), "nodes")

    # Interior beads only: element 0 and -1 of a path are the two junctions,
    # which already have their own displacement.
    bead_xyz = {}
    for k, atoms in chain_atoms.items():
        for i, idx in enumerate(atoms):
            bead_xyz[idx] = tuple(paths[k][i + 1] / s)
    write_lammps_displacement_file(bead_xyz, s, s, s,
                                   str(chem / "system_beads.displace"), "beads")

    with open(chem / "system.groups", "w") as f:
        ids = " ".join(str(i + 1) for i in sorted(node_atom.values()))
        f.write(f"group nodes id {ids}\ngroup beads subtract all nodes\n")

    return mol.GetNumAtoms(), node_atom, chain_atoms


SIM_CONFIG = {"cg": {"soft_push_steps": 20000, "ramp_steps": 20000,
                     "equil_steps": 20000}}


# Matches examples/demos/polymer/coarse_grained/basic/config.json. Not a
# tuning knob: straight chains at melt density start with bonds near 0.15
# sigma, so a cutoff of 0.85 reads every bonded neighbour as an overlap and
# pushes the chain apart. Measured at 0.85 on this network: median bond 1.53
# sigma, longest 10.9, 7227 bonds beyond 2.5, and the resolver still had
# 15679 overlaps left after 20 iterations because it was fighting the
# connectivity. At 0.01 it only separates exact coincidences and leaves the
# expansion to stage 1's soft push, which is what that stage is calibrated
# for.
CG_OVERLAP_CUTOFF = 0.01


def conform_and_script(root, graph, geo, overlap_cutoff=CG_OVERLAP_CUTOFF):
    study = root.name
    cm = ConformationManager(str(root.parent), study)
    conformed, roles = cm.apply_displacements(
        "system.data", lattice_box=tuple(geo["box"]), periodicity=(1, 1, 1))
    noisy = cm.apply_noise(conformed, magnitude=1e-4)
    cm.resolve_overlaps(noisy, roles, cutoff=overlap_cutoff, max_iters=10)

    gen = LammpsInputGenerator(str(root.parent), study, config=SIM_CONFIG)
    gen.write_serial_soft_minimization(settings_file="system.in.settings",
                                       model_type="cg")
    gen.write_parallel_production(settings_file="system.in.settings",
                                  model_type="cg")
    return root / "04_Simulation"


def run_md(sim_dir, stages=3):
    scripts = ["minimize_1_serial.in", "minimize_2_parallel.in",
               "minimize_3_parallel.in"][:stages]
    runner = SimulationRunner(sim_dir=sim_dir, executable="lmp",
                              n_procs=1, use_mpi=False)
    runner.run_sequence(scripts)
    return [sim_dir / f"log.{Path(s).stem.replace('minimize_', 'stage')}.lammps"
            for s in scripts]


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------

def step1(args):
    """The network itself. No entanglements: this is the control."""
    graph = build_network()
    geo = geometry(graph, density=args.density, bond=args.bond)

    # Step 1 has no braids, so a path is its chord and the longest path is
    # the longest chord. From step 2 on the same call re-pins the scale to
    # whatever the braids made longest, with no change here.
    if args.density is None:
        unit = linear_paths_unit(geo)
        geo = geometry(graph, bond=args.bond,
                       scale=scale_for_longest(unit, DP, args.bond))
    paths = place_beads({k: np.stack(c) for k, c in geo["chords"].items()},
                        ends=geo["ends"])

    tag = "melt" if args.density is not None else f"b{int(round(args.bond*100))}"
    root = OUT / f"step1_{tag}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)

    L = geo["L"]
    deg = sorted(dict(graph.degree()).values())
    hist = {d: deg.count(d) for d in sorted(set(deg))}
    dps = np.array([len(p) - 2 for p in paths.values()])
    built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print()
    print(f"  lattice        MIX {'x'.join(str(d) for d in LATTICE['dims'])}, "
          f"mix {LATTICE['mix']}")
    print(f"  junctions      {graph.number_of_nodes()}   functionality {hist}")
    print(f"  chains         {graph.number_of_edges()}   beads {n_atoms}")
    print(f"  scale set by   "
          f"{'density' if args.density is not None else 'bond length'}")
    print(f"  box            {L[0]:.1f} sigma cube "
          f"(junction spacing {geo['scale']:.1f})")
    print(f"  density        {geo['density']:.4f}")
    print(f"  DP             {int(np.median(dps))} for every chain")
    print(f"  bond as built  longest path {built.max():.3f}  "
          f"median {np.median(built):.3f}  shortest {built.min():.3f} sigma")
    print(f"  FENE headroom  longest bond {built.max():.3f} of 1.5 limit  "
          f"({'ok' if built.max() < 1.5 else 'OVER'})")
    print()

    sim_dir = conform_and_script(root, graph, geo)
    (root / "network.json").write_text(json.dumps(
        {"junctions": graph.number_of_nodes(), "chains": graph.number_of_edges(),
         "dp": int(np.median(dps)), "beads": n_atoms,
         "scale": geo["scale"], "density": geo["density"],
         "target_bond": args.bond,
         "built_bond_median": float(np.median(built)),
         "built_bond_max": float(built.max()),
         "box_sigma": L.tolist(), "functionality": hist}, indent=2))

    if args.run_md:
        print("--- LAMMPS ---")
        run_md(sim_dir, args.stages)
        report_bonds(root)
    else:
        print(f"  scripts in {sim_dir}  (add --run-md to run them)")
    return 0


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def read_data(path):
    """Box, coordinates and bond list from a LAMMPS data file."""
    L, xyz, bonds, sec = {}, {}, [], None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s.endswith(("xlo xhi", "ylo yhi", "zlo zhi")):
            p = s.split()
            L[p[-1][0]] = float(p[1]) - float(p[0])
        if s and s.split()[0] in ("Atoms", "Bonds", "Masses", "Velocities",
                                  "Angles"):
            sec = s.split()[0]
            continue
        if not s or s.startswith("#"):
            continue
        p = s.split()
        if sec == "Atoms" and len(p) >= 7:
            xyz[int(p[0])] = np.array([float(p[4]), float(p[5]), float(p[6])])
        elif sec == "Bonds" and len(p) == 4:
            bonds.append((int(p[2]), int(p[3])))
    return np.array([L["x"], L["y"], L["z"]]), xyz, bonds


def bond_lengths(path):
    box, xyz, bonds = read_data(path)
    out = np.empty(len(bonds))
    for i, (u, v) in enumerate(bonds):
        d = xyz[v] - xyz[u]
        out[i] = np.linalg.norm(d - box * np.round(d / box))
    return out


def report_bonds(root):
    """Bond length through the protocol. The one number that says whether
    the chains were built at the length the force field wants, or whether
    stage 2 had to stretch them there."""
    stages = [("built", "03_Conformation/system_conformed.data"),
              ("stage 1 soft", "04_Simulation/system_after_soft.data"),
              ("stage 2 ramp", "04_Simulation/system_ramped.data"),
              ("stage 3 equil", "04_Simulation/system_equilibrated.data")]
    print()
    print(f"  {'':14s} {'median':>7} {'max':>7} {'>1.5sig':>8}")
    for label, rel in stages:
        p = root / rel
        if not p.exists():
            print(f"  {label:14s}   not written")
            continue
        d = bond_lengths(p)
        print(f"  {label:14s} {np.median(d):7.3f} {d.max():7.3f} "
              f"{int((d > 1.5).sum()):7d}")


# ---------------------------------------------------------------------------
# Step 2: one pair, entangled once
# ---------------------------------------------------------------------------

def neighbour_shells(geo, max_units=0.8, tol=0.02):
    """Chord pairs grouped into the lattice's discrete separation bands.

    On a lattice the closest-approach distance between two strands takes a
    handful of values, not a continuum -- 0.20, 0.35, 0.41, 0.50, 0.61,
    0.71 lattice units for this mix. Those bands are what "first neighbour,
    second neighbour" actually means here, so shells are read off the
    geometry rather than assumed.
    """
    s, ch = geo["scale"], geo["chords"]
    ids = sorted(ch)
    found = []
    for i, ka in enumerate(ids):
        a0, a1 = ch[ka]
        mid_a = 0.5 * (a0 + a1)
        for kb in ids[i + 1:]:
            b0, b1 = ch[kb]
            if np.linalg.norm(mid_a - 0.5 * (b0 + b1)) > 2.5 * s:
                continue
            t, _ = closest_approach(a0, a1, b0, b1)
            gap, _ = gap_at(a0, a1, b0, b1, t)
            u = gap / s
            if 1e-6 < u <= max_units:
                found.append((u, ka, kb))
    found.sort()

    bands, cur = [], []
    for rec in found:
        if cur and rec[0] - cur[0][0] > tol:
            bands.append(cur)
            cur = []
        cur.append(rec)
    if cur:
        bands.append(cur)
    return bands


def build_with_braids(graph, requests, bond=BOND, dp=DP, rounds=4):
    """Paths for every chain, with the scale pinned to the longest of them.

    Braid size follows the gap and the gap follows the scale, so pinning the
    scale to a braided path is a fixed point rather than a formula. It is a
    contraction -- the braid is nearly scale-free in lattice units -- so a
    few rounds settle it.
    """
    geo = geometry(graph, bond=bond)
    alloc = None
    for _ in range(rounds):
        alloc = allocate_contacts(requests, geo["chords"], BraidShape())
        raw = {k: compose_chain_path(k, alloc, geo["chords"], 400)
               for k in geo["chords"]}
        unit = {k: p / geo["scale"] for k, p in raw.items()}
        new_scale = scale_for_longest(unit, dp, bond)
        if abs(new_scale - geo["scale"]) < 1e-6 * geo["scale"]:
            break
        geo = geometry(graph, bond=bond, scale=new_scale)

    raw = {k: compose_chain_path(k, alloc, geo["chords"], 4000)
           for k in geo["chords"]}
    return geo, alloc, place_beads(raw, dp, geo["ends"])


def unwrap_chain(ids_seq, xyz, box):
    """Bead path of one chain, unwrapped so a boundary crossing stays whole."""
    out = [xyz[ids_seq[0]]]
    for aid in ids_seq[1:]:
        d = xyz[aid] - out[-1]
        out.append(out[-1] + (d - box * np.round(d / box)))
    return np.array(out)


def chain_ids(chain, node_atom, chain_atoms, ends):
    u, v = ends[chain]
    return ([node_atom[u] + 1]
            + [i + 1 for i in chain_atoms[chain]]
            + [node_atom[v] + 1])


def measure_pair(data_file, seq_a, seq_b, contact=None):
    """Linking number and closest approach of two chains in a data file."""
    box, xyz, _ = read_data(data_file)
    pa = unwrap_chain(seq_a, xyz, box)
    pb = unwrap_chain(seq_b, xyz, box)
    # Bring b into the image nearest a, or the linking integral sees two
    # curves that never interact.
    d = pa.mean(axis=0) - pb.mean(axis=0)
    pb = pb + box * np.round(d / box)
    return far_closed_linking(pa, pb, contact), min_separation(pa, pb)


def step2(args):
    """One pair of chains, entangled once."""
    graph = build_network()
    probe = geometry(graph, bond=args.bond)
    bands = neighbour_shells(probe)
    if not bands:
        raise SystemExit("no chord pairs close enough to consider")

    print(f"  neighbour shells found (gap in lattice units):")
    for i, band in enumerate(bands[:6], start=1):
        print(f"    shell {i}: gap {band[0][0]:.2f}  ({len(band)} pairs)")

    band = bands[min(args.shell, len(bands)) - 1]
    gap_u, ka, kb = band[0]
    print(f"\n  picked shell {args.shell}: chains {ka} and {kb}, "
          f"gap {gap_u:.2f} lattice units")

    req = [ContactRequest(ka, kb, windings=args.windings)]
    geo, alloc, paths = build_with_braids(graph, req, args.bond)

    if not alloc.accepted:
        why = alloc.rejected[0].reason if alloc.rejected else "unknown"
        raise SystemExit(f"  contact refused: {why}")
    a = alloc.accepted[0]

    built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    braid_bond = np.linalg.norm(np.diff(paths[ka], axis=0), axis=1)
    lk = far_closed_linking(paths[ka], paths[kb], a.contact)
    sep = min_separation(paths[ka], paths[kb])

    t = np.linspace(0, 1, DP + 2)[:, None]
    sa = geo["chords"][ka][0] + t * (geo["chords"][ka][1] - geo["chords"][ka][0])
    sb = geo["chords"][kb][0] + t * (geo["chords"][kb][1] - geo["chords"][kb][0])
    lk0 = far_closed_linking(sa, sb, a.contact)

    print(f"\n  box            {geo['L'][0]:.1f} sigma "
          f"(spacing {geo['scale']:.1f}), density {geo['density']:.4f}")
    print(f"  contact gap    {a.contact.gap:.2f} sigma, "
          f"windings asked {args.windings}, granted {a.windings}")
    print(f"  bond as built  longest path {built.max():.3f}  "
          f"median {np.median(built):.3f} sigma   "
          f"(in the braid: max {braid_bond.max():.3f})")
    print(f"  linking number {lk0:+.2f} without the braid  ->  "
          f"{lk:+.2f} with it")
    print(f"  separation     {sep:.2f} sigma at closest approach")

    root = OUT / f"step2_shell{args.shell}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    seq_a = chain_ids(ka, node_atom, chain_atoms, geo["ends"])
    seq_b = chain_ids(kb, node_atom, chain_atoms, geo["ends"])

    sim_dir = conform_and_script(root, graph, geo)
    (root / "pair.json").write_text(json.dumps(
        {"chain_a": ka, "chain_b": kb, "shell": args.shell,
         "gap_lattice_units": gap_u, "gap_sigma": a.contact.gap,
         "windings": a.windings, "linking_before": lk0, "linking_built": lk,
         "separation_built": sep, "scale": geo["scale"],
         "density": geo["density"], "beads": n_atoms}, indent=2))

    if args.run_md:
        print("\n--- LAMMPS ---")
        run_md(sim_dir, args.stages)
        report_bonds(root)
        print()
        print(f"  {'':14s} {'linking':>8} {'separation':>11}")
        print(f"  {'as built':14s} {lk:8.2f} {sep:11.2f}")
        for label, rel in [("stage 1 soft", "04_Simulation/system_after_soft.data"),
                           ("stage 2 ramp", "04_Simulation/system_ramped.data"),
                           ("stage 3 equil", "04_Simulation/system_equilibrated.data")]:
            p = root / rel
            if not p.exists():
                continue
            l2, s2 = measure_pair(p, seq_a, seq_b, a.contact)
            # Signed, against the value as built. Comparing magnitudes would
            # call a destroyed +1 that came back as -1 a success, and that is
            # exactly what happens here: the strands cross and re-link the
            # other way round.
            keep = "kept" if abs(l2 - lk) < 0.5 else "LOST"
            print(f"  {label:14s} {l2:8.2f} {s2:11.2f}   {keep}")
    else:
        print(f"\n  scripts in {sim_dir}  (add --run-md to run them)")
    return 0


STEPS = {1: step1, 2: step2}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", type=int, default=1, choices=sorted(STEPS))
    ap.add_argument("--shell", type=int, default=1,
                    help="which neighbour shell to entangle (step 2+)")
    ap.add_argument("--windings", type=int, default=1)
    ap.add_argument("--run-md", action="store_true")
    ap.add_argument("--stages", type=int, default=3, choices=(1, 2, 3))
    ap.add_argument("--bond", type=float, default=BOND,
                    help="bond length every chain is built at, in sigma; "
                         "the scale follows from it")
    ap.add_argument("--density", type=float, default=None,
                    help="set the scale from bead density instead, e.g. 0.85 "
                         "for the melt comparison")
    args = ap.parse_args()

    print("=" * 70)
    print(f"Step {args.step}: {STEPS[args.step].__doc__.splitlines()[0]}")
    print("=" * 70)
    return STEPS[args.step](args)


if __name__ == "__main__":
    raise SystemExit(main())
