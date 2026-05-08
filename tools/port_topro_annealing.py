"""Port the 14 topro CHARMM annealing scripts to MARTINI 3.

Each anneal_NN_*.in script gets:
  * MARTINI FF setup (lj/cut/coul/cut + dielectric=15 + cosine/squared angles +
    no kspace + no fix shake on water)
  * comm_modify cutoff 60.0
  * Local-test-friendly run lengths (1000 steps x 1 fs each, instead of topro's
    500000 x 1 fs = 0.5 ns per stage). Bump `RUN_STEPS` for production.
"""
from __future__ import annotations

from pathlib import Path

# Local-test run length per anneal stage. Topro production = 500000 (= 0.5 ns).
# 1000 steps x 1 fs = 1 ps per stage -> ~14 ps total annealing for the local
# demo. Bump back to 500000 for the supercomputer run.
RUN_STEPS = 1000

_HEADER = """units           real
timestep        1.0
boundary        p p p

atom_style      full
bond_style      harmonic
angle_style     cosine/squared
dihedral_style  charmm
improper_style  harmonic

pair_style      lj/cut/coul/cut 12.0
pair_modify     shift yes
special_bonds   lj 0.0 0.0 0.0 coul 0.0 0.0 0.0
dielectric      15.0

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
comm_modify     mode single cutoff 60.0
"""

_THERMO = """thermo          {thermo_every}
thermo_style    custom step vol temp press density etotal pe ke evdwl ebond epair eangle
"""


def _stage(stage_num: int, name: str, prev_input: str, fix_block: str, run_steps: int = RUN_STEPS) -> str:
    """One anneal_NN script. `prev_input` is read first; restart is written at end."""
    is_init = (stage_num == 0)
    if is_init:
        read_block = (
            "read_data       ../system_equilibrated.data\n"
            "include         ../protein_network.in.settings\n"
            "include         ../protein_network.in.groups\n"
        )
    else:
        read_block = (
            f"read_restart    {prev_input}\n"
            "include         ../protein_network.in.groups\n"
        )
    minimize_block = ""
    velocity_block = ""
    if is_init:
        minimize_block = (
            "min_style       cg\n"
            "minimize        1.0e-5 1.0e-7 5000 50000\n"
            "write_data      anneal_00_minimized.data\n\n"
            "reset_timestep  0\n"
            "velocity        all create 310 482931 dist gaussian loop local\n"
        )
    return f"""# anneal_{stage_num:02d}_{name}.in (MARTINI 3 port of topro CHARMM)
# Reads:  {prev_input if not is_init else "../system_equilibrated.data"}
# Writes: anneal_{stage_num:02d}.restart
# RUN_STEPS = {run_steps} (= {run_steps} fs); bump to 500000 for production (0.5 ns).

{_HEADER}
{read_block}
{_THERMO.format(thermo_every=max(100, run_steps // 10))}
{minimize_block}{fix_block}
run             {run_steps}
unfix           1

write_restart   anneal_{stage_num:02d}.restart
print           "anneal_{stage_num:02d}_{name} complete"
"""


def generate_all(out_dir: Path) -> list[Path]:
    """Generate all 14 MARTINI-ported anneal_NN scripts in `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    stages = [
        (0,  "init",   None,                  "fix             1 all nvt temp 310 310 100"),
        (1,  "settle", "anneal_00.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (2,  "settle", "anneal_01.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (3,  "pramp",  "anneal_02.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 200.0 1000"),
        (4,  "tramp",  "anneal_03.restart",   "fix             1 all npt temp 310 380 100 iso 200.0 200.0 1000"),
        (5,  "hold",   "anneal_04.restart",   "fix             1 all npt temp 380 380 100 iso 200.0 200.0 1000"),
        (6,  "tcool",  "anneal_05.restart",   "fix             1 all npt temp 380 310 100 iso 200.0 200.0 1000"),
        (7,  "pdrop",  "anneal_06.restart",   "fix             1 all npt temp 310 310 100 iso 200.0 1.0 1000"),
        (8,  "relax",  "anneal_07.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (9,  "relax",  "anneal_08.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (10, "relax",  "anneal_09.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (11, "relax",  "anneal_10.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (12, "relax",  "anneal_11.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
        (13, "relax",  "anneal_12.restart",   "fix             1 all npt temp 310 310 100 iso 1.0 1.0 1000"),
    ]
    for stage_num, name, prev, fix_block in stages:
        path = out_dir / f"anneal_{stage_num:02d}_{name}.in"
        path.write_text(_stage(stage_num, name, prev, fix_block))
        paths.append(path)
    # Final write_data after the last stage
    last = paths[-1]
    last.write_text(last.read_text() + 'write_data      ../system_annealed.data\nprint           "Annealing pipeline complete -> ../system_annealed.data"\n')
    return paths


if __name__ == "__main__":
    import sys
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("annealing")
    paths = generate_all(out)
    for p in paths:
        print(f"  wrote {p}")
