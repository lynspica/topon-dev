"""Post-run summary for a `topon generate` output directory.

Used by `topon inspect <run_dir>`. Parses each stage's outputs and
prints a one-screen status report: atom counts, box, what artifacts
landed, what the next LAMMPS commands are.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


_DISPLACE_KINDS = (
    "system_nodes.displace",
    "system_backbone.displace",
    "system_beads.displace",     # legacy pre-Option-C
    "system_grafts.displace",
    "system_pendant.displace",
    "system_hydrogens.displace",
)


@dataclass
class StageReport:
    name: str
    present: bool
    files: list[str] = field(default_factory=list)
    summary: str = ""


@dataclass
class RunSummary:
    root: Path
    chemistry: StageReport
    conformation: StageReport
    simulation: StageReport
    atom_count: Optional[int] = None
    n_atom_types: Optional[int] = None
    box: Optional[tuple[float, float, float]] = None


def _parse_system_data(path: Path) -> dict:
    """Extract atom_count, n_atom_types, box from a LAMMPS data file header."""
    info: dict = {"atom_count": None, "n_atom_types": None, "box": None}
    if not path.exists():
        return info
    box_lo = [None, None, None]
    box_hi = [None, None, None]
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if i > 80:
                break
            s = line.strip()
            m = re.match(r"^(\d+)\s+atoms\s*$", s)
            if m:
                info["atom_count"] = int(m.group(1))
                continue
            m = re.match(r"^(\d+)\s+atom\s+types\s*$", s)
            if m:
                info["n_atom_types"] = int(m.group(1))
                continue
            m = re.match(r"^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([xyz])lo\s+\3hi\s*$", s)
            if m:
                idx = {"x": 0, "y": 1, "z": 2}[m.group(3)]
                box_lo[idx] = float(m.group(1))
                box_hi[idx] = float(m.group(2))
    if all(lo is not None and hi is not None for lo, hi in zip(box_lo, box_hi)):
        info["box"] = (
            box_hi[0] - box_lo[0],
            box_hi[1] - box_lo[1],
            box_hi[2] - box_lo[2],
        )
    return info


def _chemistry_summary(d: Path) -> StageReport:
    files = sorted(p.name for p in d.iterdir() if p.is_file()) if d.exists() else []
    if not files:
        return StageReport("02_Chemistry", present=False)
    data = _parse_system_data(d / "system.data")
    displace_present = [n for n in _DISPLACE_KINDS if n in files]
    bits = []
    if data["atom_count"]:
        bits.append(f"{data['atom_count']} atoms")
    if data["n_atom_types"]:
        bits.append(f"{data['n_atom_types']} atom types")
    if data["box"]:
        lx, ly, lz = data["box"]
        bits.append(f"box {lx:.1f} x {ly:.1f} x {lz:.1f} A")
    bits.append(f"{len(displace_present)} displace file(s): "
                f"{', '.join(d.replace('system_', '').replace('.displace', '') for d in displace_present)}")
    return StageReport(
        name="02_Chemistry",
        present=True,
        files=files,
        summary="; ".join(bits),
    )


def _conformation_summary(d: Path) -> StageReport:
    files = sorted(p.name for p in d.iterdir() if p.is_file()) if d.exists() else []
    if not files:
        return StageReport("03_Conformation", present=False)
    bits = []
    for name, label in (
        ("system_conformed.data", "conformed"),
        ("system_relaxed.data", "relaxed"),
    ):
        if name in files:
            data = _parse_system_data(d / name)
            if data["atom_count"]:
                bits.append(f"{label} ({data['atom_count']} atoms)")
            else:
                bits.append(label)
    return StageReport(
        name="03_Conformation",
        present=True,
        files=files,
        summary="; ".join(bits) if bits else f"{len(files)} file(s)",
    )


def _simulation_summary(d: Path) -> StageReport:
    files = sorted(p.name for p in d.iterdir() if p.is_file()) if d.exists() else []
    if not files:
        return StageReport("04_Simulation", present=False)
    scripts = [f for f in files if f.endswith(".in")]
    logs = [f for f in files if f.endswith(".lammps") or f.startswith("log.")]
    data_files = [f for f in files if f.endswith(".data")]
    bits = [
        f"{len(scripts)} LAMMPS script(s)",
        f"{len(logs)} log(s)",
        f"{len(data_files)} stage output data file(s)",
    ]
    if "system_equilibrated.data" in files:
        bits.append("equilibrated.data present")
    elif "system_minimized_final.data" in files:
        bits.append("minimized_final.data present (stage 3 minimize done)")
    elif "system_ramped.data" in files:
        bits.append("ramped.data present (stage 2 done)")
    elif "system_after_soft.data" in files:
        bits.append("after_soft.data present (stage 1 done)")
    return StageReport(
        name="04_Simulation",
        present=True,
        files=files,
        summary="; ".join(bits),
    )


def summarise(run_dir: Path) -> RunSummary:
    run_dir = Path(run_dir)

    # Layout A: study root containing 02_Chemistry/ + 03_Conformation/ + 04_Simulation/
    # Layout B: parent of A — find the nested study dir
    # Layout C: flat folder (expected_output/-style): every artifact at top level
    if (run_dir / "02_Chemistry").exists():
        root = run_dir
        flat = False
    else:
        nested = list(run_dir.glob("*/02_Chemistry"))
        if nested:
            root = nested[0].parent
            flat = False
        else:
            root = run_dir
            flat = True

    if not flat:
        chem = _chemistry_summary(root / "02_Chemistry")
        conf = _conformation_summary(root / "03_Conformation")
        sim = _simulation_summary(root / "04_Simulation")
        info = _parse_system_data(root / "02_Chemistry" / "system.data")
    else:
        # Flat: split files by name pattern into the three stage buckets.
        chem = _chemistry_summary(root)
        chem.name = "02_Chemistry (flat)"
        # No separate conformation/simulation buckets; report what's there.
        files = sorted(p.name for p in root.iterdir() if p.is_file())
        scripts = [f for f in files if f.endswith(".in")]
        logs = [f for f in files if f.endswith(".lammps") or f.startswith("log.")]
        stage_data = [f for f in files if f.endswith(".data")
                      and f not in ("system.data", "system_conformed.data",
                                    "system_relaxed.data")]
        sim_bits = []
        if scripts:
            sim_bits.append(f"{len(scripts)} LAMMPS script(s)")
        if logs:
            sim_bits.append(f"{len(logs)} log(s)")
        if stage_data:
            sim_bits.append(f"{len(stage_data)} stage data file(s)")
        if "system_equilibrated.data" in files:
            sim_bits.append("equilibrated.data present")
        elif "system_minimized_final.data" in files:
            sim_bits.append("minimized_final.data present")
        elif "system_ramped.data" in files:
            sim_bits.append("ramped.data present (stage 2 ok)")
        elif "system_after_soft.data" in files:
            sim_bits.append("after_soft.data present (stage 1 ok)")
        sim = StageReport("04_Simulation (flat)", present=bool(scripts or stage_data),
                          files=scripts + logs + stage_data,
                          summary="; ".join(sim_bits) if sim_bits else "(no LAMMPS files)")
        conf_present = any(
            n in files for n in ("system_conformed.data", "system_relaxed.data")
        )
        conf = StageReport(
            "03_Conformation (flat)",
            present=conf_present,
            files=[n for n in ("system_conformed.data", "system_relaxed.data") if n in files],
            summary="conformed+relaxed present" if conf_present else "(none)",
        )
        info = _parse_system_data(root / "system.data")

    return RunSummary(
        root=root,
        chemistry=chem,
        conformation=conf,
        simulation=sim,
        atom_count=info.get("atom_count"),
        n_atom_types=info.get("n_atom_types"),
        box=info.get("box"),
    )


def format_summary(s: RunSummary) -> str:
    """Render a RunSummary as multi-line text."""
    lines = [
        f"topon inspect  {s.root}",
        "",
        f"  atoms       : {s.atom_count if s.atom_count is not None else '?'}",
        f"  atom types  : {s.n_atom_types if s.n_atom_types is not None else '?'}",
    ]
    if s.box:
        lx, ly, lz = s.box
        lines.append(f"  box (A)     : {lx:.2f} x {ly:.2f} x {lz:.2f}")
    lines.extend([
        "",
        "  Stage status:",
    ])
    for stage in (s.chemistry, s.conformation, s.simulation):
        if stage.present:
            lines.append(f"    [ok]   {stage.name:<18s}  {stage.summary}")
        else:
            lines.append(f"    [miss] {stage.name:<18s}  (no files)")

    # Suggest the next LAMMPS command (try nested first, then flat root)
    candidates = [s.root / "04_Simulation", s.root]
    for sim_dir in candidates:
        if not sim_dir.exists():
            continue
        scripts = sorted(p.name for p in sim_dir.iterdir() if p.suffix == ".in")
        if scripts:
            lines.extend([
                "",
                "  To run LAMMPS:",
                f"    cd {sim_dir}",
            ])
            for sc in scripts:
                lines.append(f"    lmp -in {sc}")
            break
    return "\n".join(lines)
