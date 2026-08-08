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
import hashlib
import json
from collections import defaultdict
import random
import shutil
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
from topon.conformation.entanglement.waypoints import (  # noqa: E402
    Site,
    entangled_group,
    entangled_pair,
    meander_to_length,
    resample_path,
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
# Functionality up to 8, mean near 4. Twelve is reachable on this mix but it
# makes every other problem worse: a junction with twelve chains cannot seat
# their first beads a sigma apart, the mix puts some of those chains on the
# same ray so they never separate, and a braid between any two of them
# usually has a third lying in its volume.
#
# The distribution is spelled out because most mean-4 targets are not
# reachable by sculpting this lattice: of four shapes tried, three returned
# no graph. This one lands on 354 chains, mean f 4.16, max 8.
LATTICE = dict(lattice="MIX", dims=(4, 4, 4), mix={"SC": 0.2, "BCC": 0.4, "FCC": 0.4},
               cutoff=1.0, max_func=8, seed=42,
               degree_dist="0:0,1:0,2:15,3:40,4:55,5:35,6:15,7:7,8:3")
# Beads per chain between its two junctions.
#
# 80 rather than 40 because several entanglements on one pair need the two
# chains to run alongside each other for long enough to hold them. At DP 40
# the chord came out 23.5 sigma, the stretch over which the pair stayed
# close was 14.2, and one site costs 7.8 -- room for one, whatever else was
# tuned. The limit was the length of the chain, not the braid.
DP = 80
DENSITY = 0.85          # melt density, for the --density comparison run

# Bond length every chain is built at, and the number that fixes the scale:
# the longest chain, drawn straight with DP beads, gets bonds this long.
#
# The scale follows from the geometry, not from a chosen density; the density
# is reported instead. 0.90 sits just under the Kremer-Grest equilibrium of
# ~0.97, so the protocol has nothing to stretch. 0.95 rather than 0.90 buys
# a little more room for the braid to spend on clearance.
BOND = 0.95

# Minimum separation asked of the first beads of chains sharing a junction.
# Set to 0 to build without the shell, which is how the overlap counts in
# junction_shell's docstring were measured.
SHELL_SPACING = 1.0

# How much contour each chain folds into its chord: (DP+1)*BOND over the
# typical chord. This is the knob that decides whether a designed
# entanglement is distinguishable from the ones a coiled chain makes by
# accident. See geometry() for the measurement; 1.8 sits safely under the
# 2.5 where control is lost.
COIL = 1.8

# Lattice scale, in sigma per lattice unit, that BraidShape's default lengths
# were calibrated against. The braid is scaled by scale/REF_SCALE so it keeps
# the same proportions on any lattice, which is what makes the granted
# winding count independent of the box the scale search happens to be at.
REF_SCALE = 25.0

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


def build_network(spec=LATTICE, cache=True):
    """The shared network, cached to disk after the first search.

    Every step in the series has to get the *same* network or comparisons
    between them mean nothing. Seeding alone is not enough: sculpting to a
    mean-4 target on this mix is a search that succeeds on a different trial
    each run and sometimes runs out of time, so the graph is written once
    and reloaded. Delete the cache file to search again.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    # Everything that changes the network goes in the key. The mix fractions
    # and the degree target were missing, so two different mixes shared a
    # cache entry and the second silently got the first one's graph.
    mix = "-".join(f"{k}{spec['mix'].get(k, 0):.2f}"
                   for k in sorted(spec.get("mix") or {}))
    dd = hashlib.sha1(str(spec["degree_dist"]).encode()).hexdigest()[:6]
    key = (f"{spec['lattice']}_{'x'.join(str(d) for d in spec['dims'])}"
           f"_f{spec['max_func']}_s{spec['seed']}"
           f"{('_' + mix) if mix else ''}_{dd}")
    path = OUT / f"network_{key}.gpickle"

    if cache and path.exists():
        import pickle
        with open(path, "rb") as f:
            return pickle.load(f)

    random.seed(spec["seed"])
    np.random.seed(spec["seed"])
    gen = PythonTopologyGenerator(_Cfg(spec))
    graphs = gen.generate(trials=20000, max_saves=1, time_limit=900)
    if not graphs:
        raise SystemExit(
            f"sculpting produced no {spec['lattice']} network for "
            f"{spec['degree_dist']}. Not every degree distribution is "
            f"reachable on this lattice; try a broader one.")
    if cache:
        import pickle
        with open(path, "wb") as f:
            pickle.dump(graphs[0], f)
    return graphs[0]


def geometry(graph, dp=DP, density=None, bond=BOND, scale=None,
             coil=None):
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
        if coil is not None:
            # Set the scale from how much the chains have to fold up. Each
            # chain carries (dp+1)*bond of contour whatever its chord, so
            # the ratio of the two is what says whether a designed
            # entanglement can be told apart from the coil's own accidental
            # crossings. Measured on one pair, asked one, two and three
            # sites:
            #
            #   coil   1 site   2 sites  3 sites
            #   7.8x    -1.15    +1.97    -1.17
            #   5.4x    -1.45    -0.23    +0.45
            #   3.6x    +0.71    +2.43    -0.58
            #   2.5x    +1.10    +2.13    +3.09
            #   1.7x    -0.96    -2.01    -3.05
            #
            # Below about 2.5 the count is what was asked for. Above it the
            # coil makes more crossings than the design does and nothing is
            # controllable. Density is then reported, not chosen -- which is
            # the right way round, since the density here only has to let the
            # simulation run, and the topology is the thing being designed.
            c_typ = float(np.median([
                np.linalg.norm((raw[v] - raw[u])
                               - box * np.round((raw[v] - raw[u]) / box))
                for u, v in graph.edges()]))
            scale = bond * (dp + 1) / (coil * c_typ)
        elif density is not None:
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
                bond=c_max * scale / (dp + 1),
                coil=(dp + 1) * bond / (c_max * scale))


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


def longest_built_bond(paths):
    return max(float(np.linalg.norm(np.diff(p, axis=0), axis=1).max())
               for p in paths.values())


def pin_scale(graph, make_paths, bond=BOND, dp=DP, rounds=6, tol=1e-3):
    """Settle the scale so the longest bond actually built comes out at ``bond``.

    ``make_paths(geo)`` returns the finished bead paths for a geometry,
    junction shell and braids included.

    A formula is not enough here. The scale was being pinned from the chords
    alone, then the shell displaced beads afterwards and stretched the very
    bonds the rule was meant to bound -- built at a nominal 0.90 the longest
    bond came out at 1.210. Both the shell and the braid add length after the
    fact, so the rule has to be applied to what is built rather than to what
    was planned. Scaling the lattice scales those features with it, so this
    converges in a few rounds.
    """
    geo = geometry(graph, bond=bond, dp=dp)
    for _ in range(rounds):
        paths = make_paths(geo)
        built = longest_built_bond(paths)
        if abs(built - bond) <= tol * bond:
            break
        geo = geometry(graph, bond=bond, dp=dp,
                       scale=geo["scale"] * (bond / built))
    return geo, make_paths(geo)


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

# "attractive" is lj/cut 2.5, the pipeline default. "repulsive" is lj/cut
# 1.122462, the WCA truncation: purely repulsive, no cohesion. Which one is
# right depends on the density the network is built at. An attractive system
# held at rho 0.014 is far below where LJ holds together, so it contracts and
# carries the chains with it -- measured, a prescribed pair ended stage 2
# 16.73 sigma apart. WCA removes that driver.
PAIR_STYLE = "attractive"

# An alternative three-stage set that runs WCA throughout instead of the soft
# push. See tests/workflows/lammps_hardcore/README.md. Selected with
# --protocol hardcore; the generated scripts remain the default and are not
# modified.
HARDCORE_DIR = Path(__file__).resolve().parent / "lammps_hardcore"


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


def conform_and_script(root, graph, geo, overlap_cutoff=CG_OVERLAP_CUTOFF,
                       pair_style=PAIR_STYLE, protocol="generated"):
    study = root.name
    # Seed before the conformation stage. apply_noise perturbs every atom
    # with np.random and nothing else here sets the seed once the network
    # comes from cache, so two runs of an identical configuration diverged:
    # the same build gave 3/3 one time and 2/4 the next, and the difference
    # was the noise, not the design.
    np.random.seed(LATTICE["seed"])
    random.seed(LATTICE["seed"])
    cm = ConformationManager(str(root.parent), study)
    conformed, roles = cm.apply_displacements(
        "system.data", lattice_box=tuple(geo["box"]), periodicity=(1, 1, 1))
    noisy = cm.apply_noise(conformed, magnitude=1e-4)
    cm.resolve_overlaps(noisy, roles, cutoff=overlap_cutoff, max_iters=10)

    cfg = dict(SIM_CONFIG, pair_style=pair_style)
    gen = LammpsInputGenerator(str(root.parent), study, config=cfg)
    gen.write_serial_soft_minimization(settings_file="system.in.settings",
                                       model_type="cg")
    gen.write_parallel_production(settings_file="system.in.settings",
                                  model_type="cg")

    sim = root / "04_Simulation"
    if protocol == "hardcore":
        # Overwrite the generated scripts in this run's directory only. The
        # generator itself is untouched, so every other caller keeps the
        # scripts it has always had.
        for src in sorted(HARDCORE_DIR.glob("*.in")):
            shutil.copyfile(src, sim / src.name)
        print(f"  using the hard-core protocol from {HARDCORE_DIR.name}/")
    return sim


def run_md(sim_dir, stages=3):
    # Clear the stage outputs first. They are named for the stage, not the
    # run, so a shorter run leaves the previous run's later stages sitting
    # there and every later measurement reads them as if they were fresh --
    # a --stages 1 run reported stage 2 and stage 3 results it never
    # produced.
    for stale in ("min_stage_A.data", "min_stage_B.data",
                  "system_after_soft.data", "system_ramped.data",
                  "system_equilibrated.data"):
        f = Path(sim_dir) / stale
        if f.exists():
            f.unlink()

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

    # Step 1 has no braids, so a path is just its chord -- but the junction
    # shell still lengthens the bonds next to every crosslink, so the scale
    # is settled against the bonds actually built rather than against the
    # chords.
    def straight(g):
        return place_beads({k: np.stack(c) for k, c in g["chords"].items()},
                           ends=g["ends"])

    if args.density is None:
        geo, paths = pin_scale(graph, straight, args.bond)
    else:
        paths = straight(geo)

    tag = "melt" if args.density is not None else f"b{int(round(args.bond*100))}"
    tag += "" if args.pair_style == "attractive" else "_wca"
    tag += "" if args.protocol == "generated" else "_hard"
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

    sim_dir = conform_and_script(root, graph, geo,
                                 pair_style=args.pair_style,
                                 protocol=args.protocol)
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

def separation_bands(geo, max_units=1.3, tol=0.02):
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


def build_with_braids(graph, requests, bond=BOND, dp=DP, rounds=4,
                      density=None, verify=True, tolerance=0.35):
    """Paths for every chain, with the scale pinned to the longest of them.

    Braid size follows the gap and the gap follows the scale, so pinning the
    scale to a braided path is a fixed point rather than a formula. It is a
    contraction -- the braid is nearly scale-free in lattice units -- so a
    few rounds settle it.
    """
    if density is not None:
        # Density fixes the scale outright, so there is nothing to re-pin:
        # the chains carry far more contour than their chord needs and the
        # braid is paid for out of that slack rather than by stretching.
        geo = geometry(graph, density=density)
        alloc = allocate_contacts(requests, geo["chords"], BraidShape(),
                                  verify_windings=verify)
        raw = {k: compose_chain_path(k, alloc, geo["chords"], 4000)
               for k in geo["chords"]}
        return geo, alloc, place_beads(raw, dp, geo["ends"])

    holder = {}

    def braided(g):
        # Scale the braid with the lattice. BraidShape's defaults are
        # absolute lengths; left absolute they make the granted winding
        # count depend on the scale, and pin_scale changes the scale, so the
        # iteration chases itself. REF_SCALE is the lattice those defaults
        # were calibrated on.
        shape = BraidShape().scaled(g["scale"] / REF_SCALE)
        holder["alloc"] = allocate_contacts(requests, g["chords"], shape,
                                            verify_windings=verify,
                                            window_tolerance=tolerance)
        raw = {k: compose_chain_path(k, holder["alloc"], g["chords"], 4000)
               for k in g["chords"]}
        return place_beads(raw, dp, g["ends"])

    geo, paths = pin_scale(graph, braided, bond, dp, rounds=rounds)
    return geo, holder["alloc"], paths


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


def write_z1(path, chains, box):
    """Write chains in Z1 format: count, box, beads per chain, coordinates.

    Coordinates must be unwrapped -- Z1+ measures the shortest path a chain
    can take without crossing another, and a chain folded at the boundary
    is not the chain it is measuring.
    """
    lines = [f"{len(chains)}",
             f"{box[0]:.6f} {box[1]:.6f} {box[2]:.6f}",
             " ".join(str(len(c)) for c in chains)]
    for c in chains:
        lines += [f"{q[0]:.6f} {q[1]:.6f} {q[2]:.6f}" for q in c]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def z1_export(data_file, sequences, out_path):
    """Pull named chains out of a data file and write them for Z1+.

    Taking only the prescribed pair is what makes the answer unambiguous.
    Z1+ measures entanglements between whatever chains it is given, so with
    the rest of the network removed anything it reports is between these
    two and nothing else. Crosslinks go too: Z1+ analyses linear chains, and
    a junction would join several into one branched object it cannot read.
    """
    box, xyz, _ = read_data(data_file)
    chains = [unwrap_chain(seq, xyz, box) for seq in sequences]
    # Put every chain in the image nearest the first, or two strands that
    # are wound together get reported as far apart and unentangled.
    ref = chains[0].mean(axis=0)
    chains = [c + box * np.round((ref - c.mean(axis=0)) / box) for c in chains]
    return write_z1(out_path, chains, box)


def step2(args):
    """One pair of chains, entangled once."""
    graph = build_network()
    probe = geometry(graph, density=args.density, bond=args.bond)
    bands = separation_bands(probe)
    if not bands:
        raise SystemExit("no chord pairs close enough to consider")

    print(f"  separation bands found (gap in lattice units):")
    for i, band in enumerate(bands[:6], start=1):
        print(f"    band {i}: gap {band[0][0]:.2f}  ({len(band)} pairs)")

    # Take the first pair in the shell that can actually be built. Pairs in
    # one shell are all the same distance apart, so the choice among them is
    # arbitrary -- but they are not interchangeable, since a third chain may
    # lie in the braid volume of one and not the next.
    band = bands[min(args.shell, len(bands)) - 1]
    tried = {}
    for gap_u, ka, kb in band:
        req = [ContactRequest(ka, kb, windings=args.windings)]
        geo, alloc, paths = build_with_braids(graph, req, args.bond,
                                              density=args.density)
        if alloc.accepted:
            break
        why = alloc.rejected[0].reason if alloc.rejected else "unknown"
        tried[why] = tried.get(why, 0) + 1
    else:
        print(f"\n  no pair in band {args.shell} could be built "
              f"({len(band)} tried):")
        for why, n in sorted(tried.items(), key=lambda kv: -kv[1]):
            print(f"    {n:4d}  {why}")
        return 1

    if tried:
        print(f"\n  skipped {sum(tried.values())} pairs in this band:")
        for why, n in sorted(tried.items(), key=lambda kv: -kv[1]):
            print(f"    {n:4d}  {why}")
    print(f"\n  picked band {args.shell}: chains {ka} and {kb}, "
          f"gap {gap_u:.2f} lattice units")
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

    tag = ("melt" if args.density is not None else "b90")
    tag += "" if args.pair_style == "attractive" else "_wca"
    tag += "" if args.protocol == "generated" else "_hard"
    root = OUT / f"step2_band{args.shell}_{tag}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    seq_a = chain_ids(ka, node_atom, chain_atoms, geo["ends"])
    seq_b = chain_ids(kb, node_atom, chain_atoms, geo["ends"])

    sim_dir = conform_and_script(root, graph, geo,
                                 pair_style=args.pair_style,
                                 protocol=args.protocol)
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


# ---------------------------------------------------------------------------
# Z1+ driver
# ---------------------------------------------------------------------------

def _wsl_path(p):
    """Windows path to the /mnt/<drive> form WSL uses."""
    p = Path(p).resolve()
    drive = p.drive.rstrip(":").lower()
    rest = str(p)[len(p.drive):].replace("\\", "/")
    return f"/mnt/{drive}{rest}"


def run_z1(z1_dir, runner=None):
    """Z1+ over every .Z1 file in a directory. Returns {name: [Z per chain]}.

    Z1+ is not vendored: its README asks that it not be re-distributed, and
    only a Linux binary ships for the core module, so on Windows it runs
    under WSL. Returns None when it cannot be reached, since the rest of the
    step is still worth having without it.
    """
    import subprocess

    runner = runner or (Path(__file__).resolve().parent / "run_z1.sh")
    try:
        out = subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash",
             _wsl_path(runner), _wsl_path(z1_dir)],
            capture_output=True, text=True, timeout=900)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    text = out.stdout.replace("\x00", "")
    results, name = {}, None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("===") and s.endswith("==="):
            name = s.strip("= ").strip()
        elif s.startswith("per-chain Z:") and name:
            results[name] = [int(v) for v in s.split(":", 1)[1].split()]
    return results or None


# ---------------------------------------------------------------------------
# Step 3: several windings on one pair
# ---------------------------------------------------------------------------

def step3(args):
    """Several windings on one pair: does delivered track requested."""
    graph = build_network()
    probe = geometry(graph, bond=args.bond)
    bands = separation_bands(probe)
    band = bands[min(args.shell, len(bands)) - 1]

    z1_dir = OUT / "step3_z1"
    if z1_dir.exists():
        for f in z1_dir.glob("*.Z1"):
            f.unlink()
    z1_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    pair = None
    for want in range(1, args.max_windings + 1):
        # Keep the same pair across the sweep, so the only thing changing is
        # the number asked for. The first winding chooses it.
        candidates = [(g, a, b) for g, a, b in band] if pair is None else [pair]
        for gap_u, ka, kb in candidates:
            # `want` separate single-winding sites, not one braid wound
            # `want` times. Clearance is set by how tight the helix is and
            # does not depend on the count, so turns crammed into one site
            # are bought with the safety margin while extra sites cost only
            # chord. Allocation spreads them along the shared stretch.
            reqs = [ContactRequest(ka, kb, windings=1, priority=-i)
                    for i in range(want)]
            geo, alloc, paths = build_with_braids(graph, reqs, args.bond)
            if alloc.accepted:
                pair = (gap_u, ka, kb)
                break
        else:
            rows.append((want, None, None, None, None, "no pair could be built"))
            continue

        sites = len(alloc.accepted)
        total = sum(x.windings for x in alloc.accepted)
        built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                                for p in paths.values()])
        sep = min_separation(paths[ka], paths[kb])

        root = OUT / f"step3_e{want}"
        n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
        seq_a = chain_ids(ka, node_atom, chain_atoms, geo["ends"])
        seq_b = chain_ids(kb, node_atom, chain_atoms, geo["ends"])
        sim_dir = conform_and_script(root, graph, geo,
                                     pair_style=args.pair_style,
                                     protocol=args.protocol)
        z1_export(root / "03_Conformation/system_conformed.data",
                  [seq_a, seq_b], z1_dir / f"e{want}_1built.Z1")

        if args.run_md:
            run_md(sim_dir, args.stages)
            final = root / "04_Simulation/system_equilibrated.data"
            if final.exists():
                z1_export(final, [seq_a, seq_b], z1_dir / f"e{want}_2equil.Z1")

        # Where along the chain each site landed, as a fraction of the chord.
        c0, c1 = geo["chords"][ka]
        chord = c1 - c0
        at = sorted(float((x.contact.origin - c0) @ chord / (chord @ chord))
                    for x in alloc.accepted)
        rows.append((want, total, sites, at, sep, built.max(), None))
        print(f"  asked {want}: {sites} sites carrying {total} windings at "
              f"{', '.join(f'{v:.2f}' for v in at)} along the chain; "
              f"clearance {sep:.2f}")

    print()
    print(f"  running Z1+ on {len(list(z1_dir.glob('*.Z1')))} configurations...")
    z = run_z1(z1_dir)

    print()
    print("  " + "-" * 66)
    print(f"  {'asked':>6} {'sites':>6} {'total':>6} {'Z built':>9} "
          f"{'Z equil':>9} {'clearance':>10} {'placed at':>20}")
    print("  " + "-" * 66)
    for want, total, sites, at, sep, mx, err in rows:
        if err:
            print(f"  {want:6d}   {err}")
            continue
        zb = z.get(f"e{want}_1built") if z else None
        ze = z.get(f"e{want}_2equil") if z else None
        fb = "/".join(str(v) for v in zb) if zb else "  -"
        fe = "/".join(str(v) for v in ze) if ze else "  -"
        where = ",".join(f"{v:.2f}" for v in at)
        print(f"  {want:6d} {sites:6d} {total:6d} {fb:>9} {fe:>9} "
              f"{sep:10.2f} {where:>20}")
    print("  " + "-" * 66)
    if z is None:
        print("  Z1+ unavailable; the .Z1 files are in", z1_dir.name)
    print()
    print("  Z is entanglements per chain, one number per partner. Two chains")
    print("  wound around each other e times should each report e.")
    return 0


# ---------------------------------------------------------------------------
# Step 4: entanglements placed by hand, checked in the written data file
# ---------------------------------------------------------------------------

def build_with_waypoints(graph, pair, sites, bond=BOND, dp=DP, reach=0.45,
                         fene=1.5, density=None, coil=COIL):
    """Every chain's path, with one pair wound at the sites given.

    The sites are positions along the chain, chosen by the caller. Nothing
    is searched for and nothing is refused; see
    topon.conformation.entanglement.waypoints.

    The scale is pinned to the *straight* chains, not to the longest path.
    Pinning the longest is right when every chain is a chord, and wrong once
    one of them carries detours: an entangled chain is far longer than any
    plain chord, so pinning it shrinks the box until all the others are
    crushed. Measured with three sites on one pair, scale pinned to it: the
    median chain sat at 0.146 sigma bonds, and stage 1 had to expand 354 of
    them, leaving 438 past the FENE limit and a longest bond of 2.488.

    The entangled chain then pays for its detours out of the slack it
    already has -- at melt-like scales a chain has several times more
    contour than its chord needs -- and ``reach`` is reduced if that is not
    enough to keep it under ``fene``. Bead count stays the same for every
    chain, so DP still means one thing.
    """
    ka, kb = pair

    # Every chain is drawn at the length its beads need, entangled or not.
    #
    # Drawing the plain chains straight and sizing the box so the longest of
    # them lands at `bond` is what forced everything else: it fully extends
    # the longest chord, so on a lattice whose chords are all identical --
    # SC -- every chain is extended, no chain has slack, and an entanglement
    # has nothing to detour with. Measured there: clearance 0.00, longest
    # bond 2.843, and the pair reading 1/1 as built and 0/0 after stage 1.
    #
    # Once each path carries its own length the scale stops being determined
    # by bond length at all, and density becomes a free choice again.
    geo = geometry(graph, dp=dp, density=density, bond=bond,
                   coil=coil)

    target = (dp + 1) * bond

    def plain(c0, c1):
        """A chain with no entanglement, still drawn at the right length."""
        dense = resample(np.stack([c0, c1]), max(6 * dp, 600))
        return resample_path(meander_to_length(dense, target), dp + 2)

    def build(r):
        out = {k: plain(c0, c1)
               for k, (c0, c1) in geo["chords"].items() if k not in (ka, kb)}
        a0, a1 = geo["chords"][ka]
        b0, b1 = geo["chords"][kb]
        # bond= draws each path at the length its beads need, so the
        # entangled chain comes out at the same spacing as every other one
        # whatever its detours cost.
        pa, pb, nfo = entangled_pair(a0, a1, b0, b1, sites,
                                     n_beads=dp + 2, reach=r, bond=bond)
        out[ka], out[kb] = pa, pb
        out = apply_junction_shells(out, geo["ends"], spacing=SHELL_SPACING)
        worst = max(float(np.linalg.norm(np.diff(out[c], axis=0), axis=1).max())
                    for c in (ka, kb))
        return out, nfo, worst

    paths, info, worst = build(reach)
    return geo, paths, info, info[0].get("reach", reach)


def parse_sites(at, count, turns, span):
    """Sites from an explicit list of positions, or evenly spread.

    Explicit positions are the point of the waypoint construction: placement
    is meant to be heterogeneous, and "evenly spread" is only the fallback
    for when the caller has not said otherwise.
    """
    if at:
        spec = []
        for item in at:
            if ":" in item:
                pos, t = item.split(":", 1)
                spec.append(Site(at=float(pos), turns=int(t), span=span))
            else:
                spec.append(Site(at=float(item), turns=turns, span=span))
        return sorted(spec, key=lambda s: s.at)
    n = max(1, count)
    return [Site(at=(i + 1) / (n + 1), turns=turns, span=span)
            for i in range(n)]


def step4(args):
    """Entanglements placed by hand, checked in the written data file."""
    graph = build_network()
    probe = geometry(graph, bond=args.bond)
    bands = separation_bands(probe)
    band = bands[min(args.shell, len(bands)) - 1]
    _, ka, kb = band[0]

    sites = parse_sites(args.at, args.sites, args.windings, args.span)
    asked = sum(s.turns for s in sites)
    print(f"  pair {ka}-{kb}: "
          + ", ".join(f"{s.turns} turn(s) at {s.at:.2f}" for s in sites)
          + f"  ({asked} in total)")

    geo, paths, info, reach_used = build_with_waypoints(
        graph, (ka, kb), sites, args.bond, DP, args.reach,
        density=args.density, coil=args.coil)
    straight_a = resample(np.stack(geo["chords"][ka]), DP + 2)
    straight_b = resample(np.stack(geo["chords"][kb]), DP + 2)
    base = far_closed_linking(straight_a, straight_b)
    lk = far_closed_linking(paths[ka], paths[kb]) - base
    sep = min_separation(paths[ka], paths[kb])
    built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"  box {geo['L'][0]:.1f} sigma, density {geo['density']:.4f}")
    print(f"  as built: linking {lk:+.2f}, clearance {sep:.2f}, "
          f"longest bond {built.max():.3f}")

    root = OUT / f"step4_{LATTICE['lattice']}_{len(sites)}sites"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    seq_a = chain_ids(ka, node_atom, chain_atoms, geo["ends"])
    seq_b = chain_ids(kb, node_atom, chain_atoms, geo["ends"])
    sim_dir = conform_and_script(root, graph, geo,
                                 pair_style=args.pair_style,
                                 protocol=args.protocol)

    z1_dir = OUT / "step4_z1"
    if z1_dir.exists():
        for f in z1_dir.glob("*.Z1"):
            f.unlink()
    z1_dir.mkdir(parents=True, exist_ok=True)
    z1_export(root / "03_Conformation/system_conformed.data",
              [seq_a, seq_b], z1_dir / "1built.Z1")

    if args.run_md:
        print(f"\n--- LAMMPS, {args.stages} stage(s) ---")
        run_md(sim_dir, args.stages)
        for tag, rel in (("2stage1", "04_Simulation/system_after_soft.data"),
                         ("3stage2", "04_Simulation/system_ramped.data"),
                         ("4stage3", "04_Simulation/system_equilibrated.data")):
            out_file = root / rel
            if out_file.exists():
                z1_export(out_file, [seq_a, seq_b], z1_dir / f"{tag}.Z1")
        print()
        report_bonds(root)

    print()
    z = run_z1(z1_dir)
    print("  " + "-" * 52)
    print(f"  {'':16s} {'asked':>7} {'Z per chain':>14}")
    print("  " + "-" * 52)
    for label, key in (("as built", "1built"), ("after stage 1", "2stage1"),
                       ("after stage 2", "3stage2"),
                       ("after stage 3", "4stage3")):
        got = z.get(key) if z else None
        if got is None and key != "1built" and not (z1_dir / f"{key}.Z1").exists():
            continue
        shown = "/".join(str(v) for v in got) if got else "-"
        print(f"  {label:16s} {asked:7d} {shown:>14}")
    print("  " + "-" * 52)
    if z is None:
        print(f"  Z1+ unavailable; .Z1 files are in {z1_dir.name}")
    print()
    print("  Z is entanglements per chain, measured on the written data file")
    print("  with the rest of the network removed, so it counts only these two.")
    return 0


# ---------------------------------------------------------------------------
# Step 5: a plan of several pairs, chains carrying several partners
# ---------------------------------------------------------------------------

def build_group(graph, plan, bond=BOND, dp=DP, reach=0.45, coil=COIL):
    """Every chain in the network, with ``plan``'s pairs entangled.

    ``plan`` is ``[(chain_a, chain_b, [Site, ...]), ...]``. A chain may
    appear in several entries; it gets one path carrying all of them.
    """
    geo = geometry(graph, dp=dp, bond=bond, coil=coil)
    named = {k for a, b, _ in plan for k in (a, b)}
    sub = {k: geo["chords"][k] for k in named}
    dropped = []
    ent, info = entangled_group(sub, plan, n_beads=dp + 2, reach=reach,
                                bond=bond, dropped=dropped)
    build_group.dropped = dropped

    target = (dp + 1) * bond
    paths = dict(ent)
    for k, (c0, c1) in geo["chords"].items():
        if k in named:
            continue
        dense = resample(np.stack([c0, c1]), max(6 * dp, 600))
        paths[k] = resample_path(meander_to_length(dense, target), dp + 2)

    paths = apply_junction_shells(paths, geo["ends"], spacing=SHELL_SPACING)
    return geo, paths, info


def pick_plan(geo, spec, partners=2):
    """Turn a short description into a concrete plan on this network.

    ``spec`` is one of:
      "hub"       one chain entangled with as many close partners as it has
      "chain"     A-B, B-C, C-D, a run of pairs
      "composite" A-B, A-C, D-B twice, D-E once -- a chain that is both a
                  requester and someone else's partner, which is the case
                  the whole allocation question was about

    Partners are drawn from the closest separation band. That is not a
    detail: the contour a site costs is set by how far the two chains have
    to reach for each other, so distant partners are what make a chain run
    out of slack. Measured, one chain with two partners: from the closest
    band it builds at a coil of 1.8, from bands one to three it needs 2.5,
    which is already at the limit where control is lost.
    """
    from collections import defaultdict

    bands = separation_bands(geo)
    close = bands[0] + (bands[1] if len(bands) > 1 else [])
    # Keep the gap with each neighbour: a site costs contour in proportion
    # to how far the two chains reach for each other, so which partners a
    # hub takes decides whether it can afford them at all. Measured, one
    # chain with two partners: the two closest fit at a coil of 1.8, an
    # arbitrary two of the same neighbours need more than 91 sigma against
    # the 77 the chain has.
    nbrs = defaultdict(list)
    for gap, ka, kb in close:
        nbrs[ka].append((gap, kb))
        nbrs[kb].append((gap, ka))
    for k in nbrs:
        nbrs[k].sort()

    def sites(n):
        return [Site(at=(i + 1) / (n + 1), turns=1) for i in range(n)]

    if spec == "hub":
        # Cheapest hub, not the busiest one. A site costs contour in
        # proportion to the gap it bridges, and a chain has only so much, so
        # the chain that can afford several partners is the one whose
        # nearest few are nearest -- not the one with the longest list.
        # Measured: the busiest chain's own two closest partners still
        # needed 85 sigma before any winding, against the 77 it had.
        able = [k for k in nbrs if len(nbrs[k]) >= partners]
        if not able:
            raise SystemExit(f"no chain has {partners} close partners")
        hub = min(able, key=lambda k: sum(g for g, _ in nbrs[k][:partners]))
        picks = [p for _, p in nbrs[hub][:partners]]
        n = len(picks)
        return [(hub, p, [Site(at=(i + 1) / (n + 1), turns=1)])
                for i, p in enumerate(picks)], f"hub {hub} with {picks}"

    if spec == "chain":
        start = max(nbrs, key=lambda k: len(nbrs[k]))
        run, seen = [start], {start}
        while len(run) < 4:
            nxt = [p for _, p in nbrs.get(run[-1], []) if p not in seen]
            if not nxt:
                break
            run.append(nxt[0])
            seen.add(nxt[0])
        return ([(run[i], run[i + 1], [Site(0.5, 1)])
                 for i in range(len(run) - 1)],
                "run " + "-".join(str(x) for x in run))

    # composite: A-B, A-C, D-B, D-E, all single sites.
    #
    # B is both A's partner and D's partner, which is the case the whole
    # allocation question was about: a chain that appears in someone else's
    # plan while carrying its own. D-B was originally to be asked for twice,
    # but two sites on one pair need chains that run alongside each other and
    # this lattice's close pairs cross instead -- so every pair here gets one.
    #
    # Chosen by the worst chain, not the total. A site costs contour in
    # proportion to the gap it bridges, and every chain has the same bead
    # count, so what decides whether a plan can be built is the busiest
    # chain's own bill -- here A and D, each carrying two partners. Summing
    # over all four pairs picked a quintet needing 98.8 sigma against the 77
    # available, because a cheap pair elsewhere hid an expensive one.
    best = None
    for A in nbrs:
        if len(nbrs[A]) < 2:
            continue
        (gab, B), (gac, C) = nbrs[A][0], nbrs[A][1]
        for gdb, D in nbrs[B]:
            if D in (A, C) or len(nbrs.get(D, [])) < 2:
                continue
            rest = [(gg, k) for gg, k in nbrs[D] if k not in (A, B, C)]
            if not rest:
                continue
            gde, E = rest[0]
            cost = max(gab + gac, gdb + gde)
            if best is None or cost < best[0]:
                best = (cost, A, B, C, D, E)
    if best is None:
        raise SystemExit("no quintet supports the composite shape")

    _, A, B, C, D, E = best
    # A chain with two partners meets them at different points along
    # itself. Putting both at the midpoint makes it detour to two places at
    # once, which costs as much as one very long detour and reads as one
    # entanglement: the same quintet needed 98.8 sigma that way against the
    # 77 it had.
    plan = [(A, B, [Site(0.35, 1)]),
            (A, C, [Site(0.65, 1)]),
            (D, B, [Site(0.35, 1)]),
            (D, E, [Site(0.65, 1)])]
    label = f"{A}-{B}, {A}-{C}, {D}-{B}, {D}-{E}  (B={B} serves two)"
    return plan, label


def step5(args):
    """A plan of several pairs, with chains carrying several partners."""
    graph = build_network()
    probe = geometry(graph, dp=DP, bond=args.bond, coil=args.coil)
    plan, label = pick_plan(probe, args.plan, args.partners)
    print(f"  plan ({args.plan}): {label}")
    for a, b, sites in plan:
        print(f"    {a}-{b}: " + ", ".join(
            f"{s.turns} turn(s) at {s.at:.2f}" for s in sites))

    wanted = {}
    for a, b, sites in plan:
        n = sum(s.turns for s in sites)
        wanted[a] = wanted.get(a, 0) + n
        wanted[b] = wanted.get(b, 0) + n

    try:
        geo, paths, info = build_group(graph, plan, args.bond, DP,
                                       args.reach, args.coil)
    except ValueError as exc:
        print(f"\n  {exc}")
        return 1

    print(f"\n  box {geo['L'][0]:.1f} sigma, density {geo['density']:.4f}, "
          f"coil {args.coil}, reach solved to {info[0]['reach']:.2f}")
    built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"  bonds {built.min():.3f} to {built.max():.3f}")

    # Two ways a site can be built and still not be an entanglement, both
    # worth saying out loud rather than leaving to be discovered in the Z1+
    # table: it can land off the partner's chain, or the two chains can end
    # up lying against each other instead of winding.
    for d in info:
        a, b = d["pair"]
        if d.get("off_end"):
            print(f"  warning: the {a}-{b} site at {d['at_a']:.2f} along {a} "
                  f"falls at {d['at_b_raw']:+.2f} along {b}, past its end; "
                  f"clamped to {d['at_b']:.2f}")
    for a, b, _ in plan:
        sep = min_separation(paths[a], paths[b])
        if sep < 0.5:
            print(f"  warning: {a} and {b} are built {sep:.2f} sigma apart, "
                  f"which is lying against each other rather than winding")

    root = OUT / f"step5_{args.plan}{args.partners}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    chains = sorted(wanted)
    seqs = [chain_ids(k, node_atom, chain_atoms, geo["ends"]) for k in chains]
    sim_dir = conform_and_script(root, graph, geo,
                                 pair_style=args.pair_style,
                                 protocol=args.protocol)

    z1_dir = OUT / f"step5_{args.plan}{args.partners}_z1"
    if z1_dir.exists():
        for f in z1_dir.glob("*.Z1"):
            f.unlink()
    z1_dir.mkdir(parents=True, exist_ok=True)
    # Every pair on its own as well as all of them together. Z1+ counts
    # entanglements among whatever chains it is given, so a file holding
    # three chains cannot say whether a chain's count came from the partner
    # it was asked to entangle or from the other one wandering past. One
    # file per pair answers that; the combined file gives the total.
    seq_of = dict(zip(chains, seqs))
    want_pair = {}
    for a, b, sites in plan:
        key = (min(a, b), max(a, b))
        want_pair[key] = want_pair.get(key, 0) + sum(s.turns for s in sites)
    pairs = [(a, b) for i, a in enumerate(chains) for b in chains[i + 1:]]

    stage_files = [("1built", "03_Conformation/system_conformed.data")]
    if args.run_md:
        print(f"\n--- LAMMPS, {args.stages} stage(s) ---")
        run_md(sim_dir, args.stages)
        stage_files += [
            ("2stage1", "04_Simulation/system_after_soft.data"),
            ("3stage2", "04_Simulation/system_ramped.data"),
            ("4stage3", "04_Simulation/system_equilibrated.data")]

    for tag, rel in stage_files:
        out_file = root / rel
        if not out_file.exists():
            continue
        z1_export(out_file, seqs, z1_dir / f"{tag}_all.Z1")
        for a, b in pairs:
            z1_export(out_file, [seq_of[a], seq_of[b]],
                      z1_dir / f"{tag}_p{a}_{b}.Z1")
    if args.run_md:
        print()
        report_bonds(root)

    print()
    z = run_z1(z1_dir)
    print("  " + "-" * 60)
    head = "  ".join(f"{k:>5}" for k in chains)
    print(f"  {'':16s} {head}")
    print(f"  {'asked':16s} " + "  ".join(f"{wanted[k]:5d}" for k in chains))
    print("  " + "-" * 60)
    for lbl, key in (("as built", "1built"), ("after stage 1", "2stage1"),
                     ("after stage 2", "3stage2"),
                     ("after stage 3", "4stage3")):
        if not (z1_dir / f"{key}_all.Z1").exists():
            continue
        got = z.get(f"{key}_all") if z else None
        cells = ("  ".join(f"{v:5d}" for v in got) if got
                 else "  ".join("    -" for _ in chains))
        print(f"  {lbl:16s} {cells}")
    print("  " + "-" * 60)
    print("  Z per chain with all the named chains in one file: the total,")
    print("  including anything two of them do to each other unasked.")

    print()
    print("  " + "-" * 60)
    print(f"  {'pair on its own':20s} {'asked':>6}  " + "  ".join(
        f"{lbl:>9}" for lbl in ("built", "stage 1", "stage 2", "stage 3")))
    print("  " + "-" * 60)
    for a, b in pairs:
        cells = []
        for key in ("1built", "2stage1", "3stage2", "4stage3"):
            got = z.get(f"{key}_p{a}_{b}") if z else None
            cells.append("/".join(str(v) for v in got) if got else "-")
        want = want_pair.get((a, b))
        tag = f"{a}-{b}" if want is not None else f"{a}-{b} (not asked)"
        print(f"  {tag:20s} {want or 0:6d}  "
              + "  ".join(f"{c:>9}" for c in cells))
    print("  " + "-" * 60)
    return 0


# ---------------------------------------------------------------------------
# Step 6: shell-weighted selection across the whole network
# ---------------------------------------------------------------------------

def shell_plan(geo, weights, total, max_partners=1, seed=42):
    """Pairs drawn from named separation bands in prescribed proportions.

    ``weights`` maps band number to its share, e.g. ``{1: 0.4, 2: 0.4,
    3: 0.2}``. This is the thing the whole construction was for: a network
    where first, second and third neighbours are entangled at rates you
    choose rather than at whatever rate the geometry happens to produce.

    Bands are read off the lattice rather than assumed -- on this mix the
    closest-approach distance between two strands takes discrete values, and
    those bands are what "first neighbour" means here.

    A band can run out. Band 1 holds 8 pairs against band 5's 609, so a
    share that asks for more than a band has is served short and the
    shortfall reported rather than quietly made up from elsewhere.
    """
    rng = random.Random(seed)
    bands = separation_bands(geo)
    want = {b: int(round(total * w / sum(weights.values())))
            for b, w in weights.items()}

    used = defaultdict(int)
    plan, got, short = [], defaultdict(int), {}
    for b in sorted(want):
        if b > len(bands):
            short[b] = want[b]
            continue
        pool = list(bands[b - 1])
        rng.shuffle(pool)
        for _, ka, kb in pool:
            if got[b] >= want[b]:
                break
            if used[ka] >= max_partners or used[kb] >= max_partners:
                continue
            plan.append((ka, kb, [Site(0.5, 1)]))
            used[ka] += 1
            used[kb] += 1
            got[b] += 1
        if got[b] < want[b]:
            short[b] = want[b] - got[b]
    return plan, want, dict(got), short


def step6(args):
    """Shell-weighted entanglement across the whole network."""
    graph = build_network()
    geo0 = geometry(graph, dp=DP, bond=args.bond, coil=args.coil)
    weights = {}
    for item in (args.weights or ["1:0.4", "2:0.4", "3:0.2"]):
        b, w = item.split(":")
        weights[int(b)] = float(w)

    plan, want, got, short = shell_plan(geo0, weights, args.total,
                                        args.max_partners)
    print(f"  weights {weights}, asked {args.total} entanglements")
    for b in sorted(want):
        note = f"  ({short[b]} short, the band ran out)" if b in short else ""
        print(f"    band {b} (gap {separation_bands(geo0)[b-1][0][0]:.2f} "
              f"lattice units): asked {want[b]}, placed {got.get(b, 0)}{note}")
    if not plan:
        print("  nothing to build")
        return 1
    print(f"  {len(plan)} pairs over "
          f"{len({k for a, b, _ in plan for k in (a, b)})} chains")

    try:
        geo, paths, info = build_group(graph, plan, args.bond, DP,
                                       args.reach, args.coil)
    except ValueError as exc:
        print(f"\n  {exc}")
        return 1

    built = np.concatenate([np.linalg.norm(np.diff(p, axis=0), axis=1)
                            for p in paths.values()])
    print(f"\n  box {geo['L'][0]:.1f} sigma, density {geo['density']:.4f}, "
          f"reach solved to {info[0]['reach']:.2f}")
    print(f"  bonds {built.min():.3f} to {built.max():.3f}")
    tight = [(a, b) for a, b, _ in plan
             if min_separation(paths[a], paths[b]) < 0.5]
    if tight:
        print(f"  {len(tight)} of {len(plan)} pairs are built closer than "
              f"0.5 sigma, which is lying together rather than winding")

    root = OUT / f"step6_{args.total}"
    n_atoms, node_atom, chain_atoms = write_system(graph, geo, paths, root)
    sim_dir = conform_and_script(root, graph, geo,
                                 pair_style=args.pair_style,
                                 protocol=args.protocol)

    z1_dir = OUT / f"step6_{args.total}_z1"
    if z1_dir.exists():
        for f in z1_dir.glob("*.Z1"):
            f.unlink()
    z1_dir.mkdir(parents=True, exist_ok=True)

    # Each pair on its own. At this scale that is the only honest measure:
    # one file holding every entangled chain reports a total that cannot be
    # attributed, and the question is whether each prescribed pair got what
    # it was asked for.
    seq = {k: chain_ids(k, node_atom, chain_atoms, geo["ends"])
           for k in {c for a, b, _ in plan for c in (a, b)}}

    stage_files = [("built", "03_Conformation/system_conformed.data")]
    if args.run_md:
        print(f"\n--- LAMMPS, {args.stages} stage(s) ---")
        run_md(sim_dir, args.stages)
        stage_files = [("final", {1: "04_Simulation/system_after_soft.data",
                                  2: "04_Simulation/system_ramped.data",
                                  3: "04_Simulation/system_equilibrated.data"
                                  }[args.stages])]
        print()
        report_bonds(root)

    for tag, rel in stage_files:
        out_file = root / rel
        if not out_file.exists():
            continue
        for i, (a, b, _) in enumerate(plan):
            z1_export(out_file, [seq[a], seq[b]], z1_dir / f"{tag}_{i:03d}.Z1")

    print(f"\n  measuring {len(plan)} pairs with Z1+ ...")
    z = run_z1(z1_dir)
    if z is None:
        print(f"  Z1+ unavailable; files in {z1_dir.name}")
        return 0

    tag = stage_files[0][0]
    hit = over = under = 0
    by_band = defaultdict(lambda: [0, 0])
    bands = separation_bands(geo0)
    band_of = {}
    for bi, band in enumerate(bands, 1):
        for _, ka, kb in band:
            band_of[(min(ka, kb), max(ka, kb))] = bi
    for i, (a, b, _) in enumerate(plan):
        v = z.get(f"{tag}_{i:03d}")
        if not v:
            continue
        bi = band_of.get((min(a, b), max(a, b)), 0)
        by_band[bi][1] += 1
        if all(x == 1 for x in v):
            hit += 1
            by_band[bi][0] += 1
        elif any(x > 1 for x in v):
            over += 1
        else:
            under += 1

    print()
    print("  " + "-" * 52)
    print(f"  {len(plan)} prescribed entanglements, each pair measured alone")
    print("  " + "-" * 52)
    print(f"  exactly one, as asked   {hit:5d}   {100.0*hit/max(len(plan),1):5.1f}%")
    print(f"  more than one           {over:5d}")
    print(f"  none                    {under:5d}")
    print("  " + "-" * 52)
    for bi in sorted(by_band):
        ok, tot = by_band[bi]
        print(f"  band {bi}: {ok} of {tot} exact")
    print("  " + "-" * 52)
    return 0


STEPS = {1: step1, 2: step2, 3: step3, 4: step4, 5: step5, 6: step6}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", type=int, default=1, choices=sorted(STEPS))
    ap.add_argument("--shell", type=int, default=1,
                    metavar="BAND",
                    help="which separation band to entangle (step 2+): 1 is "
                         "the closest pair of strands, 2 the next. Distinct "
                         "from the junction shell, which is how chains leave "
                         "a crosslink.")
    ap.add_argument("--windings", type=int, default=1)
    ap.add_argument("--max-windings", type=int, default=4,
                    help="step 3 sweeps 1 up to this")
    ap.add_argument("--sites", type=int, default=3,
                    help="step 4: how many sites, spread evenly, when --at "
                         "is not given")
    ap.add_argument("--at", nargs="*", default=None, metavar="POS[:TURNS]",
                    help="step 4: place sites explicitly along the chain, "
                         "e.g. --at 0.15 0.5:2 0.85. Positions are fractions "
                         "of the chain and need not be evenly spread")
    ap.add_argument("--total", type=int, default=40,
                    help="step 6: how many entanglements to place")
    ap.add_argument("--weights", nargs="*", default=None,
                    metavar="BAND:SHARE",
                    help="step 6: how to split them across separation "
                         "bands, e.g. --weights 1:0.4 2:0.4 3:0.2")
    ap.add_argument("--max-partners", type=int, default=1,
                    help="step 6: how many entanglements one chain may carry")
    ap.add_argument("--partners", type=int, default=2,
                    help="step 5: how many partners the hub carries. Each "
                         "costs contour, so more needs a larger coil")
    ap.add_argument("--plan", default="hub",
                    choices=("hub", "chain", "composite"),
                    help="step 5: which shape of plan to build")
    ap.add_argument("--coil", type=float, default=COIL,
                    help="contour over chord, the knob that sets whether a "
                         "designed entanglement can be told apart from the "
                         "coil's own crossings. Above ~2.5 it cannot")
    ap.add_argument("--span", type=float, default=None,
                    help="step 4: how much of the chain each site occupies, "
                         "as a fraction. Smaller keeps neighbouring sites "
                         "from merging")
    ap.add_argument("--lattice", default=None,
                    help="override the lattice, e.g. SC")
    ap.add_argument("--reach", type=float, default=0.45,
                    help="step 4: how far each chain swings toward its "
                         "partner, as a fraction of their gap")
    ap.add_argument("--protocol", default="generated",
                    choices=("generated", "hardcore"),
                    help="generated uses the pipeline's own scripts; "
                         "hardcore uses tests/workflows/lammps_hardcore, "
                         "which runs WCA throughout instead of the soft push")
    ap.add_argument("--pair-style", default=PAIR_STYLE,
                    choices=("attractive", "repulsive"),
                    help="attractive is lj/cut 2.5 (pipeline default); "
                         "repulsive is lj/cut 1.122462, i.e. WCA")
    ap.add_argument("--run-md", action="store_true")
    ap.add_argument("--stages", type=int, default=3, choices=(1, 2, 3))
    ap.add_argument("--bond", type=float, default=BOND,
                    help="bond length every chain is built at, in sigma; "
                         "the scale follows from it")
    ap.add_argument("--density", type=float, default=None,
                    help="set the scale from bead density instead, e.g. 0.85 "
                         "for the melt comparison")
    args = ap.parse_args()

    if args.lattice:
        LATTICE["lattice"] = args.lattice
        if args.lattice.upper() == "SC":
            # SC 4x4x4 has 64 sites against the mix's 170, so the elaborate
            # spec the mix needs is not reachable here -- every explicit
            # mean-4 shape tried returned no graph. Capping functionality at
            # 4 and forbidding dangling ends is enough on its own and gives
            # {2:11, 3:22, 4:31}, mean 3.31.
            #
            # The geometry is the reason to try SC at all: every chord is
            # one lattice unit along an axis, so neighbouring strands are
            # parallel and evenly spaced, with none of the mix's collinear
            # pairs or wildly varying gaps.
            LATTICE["max_func"] = 4
            LATTICE["degree_dist"] = "0:0,1:0"

    print("=" * 70)
    print(f"Step {args.step}: {STEPS[args.step].__doc__.splitlines()[0]}")
    print("=" * 70)
    return STEPS[args.step](args)


if __name__ == "__main__":
    raise SystemExit(main())
