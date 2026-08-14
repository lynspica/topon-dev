"""Characteristic-ratio figure in the fig2d house style.

    python tests/workflows/char12_fig.py

C(s) only -- the internal-distance panels are dropped. Two panels, DP 40 and
DP 80, six lattices each, styled after
stage_2_swelling_design/figures/fig2d_swelling_fA_diblock.png: heavy lines,
open markers, big type, light grid, full box frame, annotations in the plot
rather than a caption.

Numbers come from the user's characteristic_ratio.py, imported from its own
location.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = Path("E:/PhD/topology_datasets")
OUTDIR = ROOT / "tests/output/char12"

sys.path.insert(0, str(REF))
from characteristic_ratio import char_ratio  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CASES = [
    ("SC",        "SC",           "#2C3E50"),
    ("BCC",       "BCC",          "#B03A2E"),
    ("FCC",       "FCC",          "#3A76B0"),
    ("MIX801010", "MIX 80/10/10", "#0e7c6b"),
    ("MIXeven",   "MIX 34/33/33", "#8E6C3A"),
    ("MIX101080", "MIX 10/10/80", "#7B4FA0"),
]

plt.rcParams.update({
    "font.size": 21,
    "axes.linewidth": 2.2,
    "axes.labelsize": 25,
    "axes.titlesize": 24,
    "xtick.labelsize": 21,
    "ytick.labelsize": 21,
    "xtick.major.width": 2.2,
    "ytick.major.width": 2.2,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "legend.fontsize": 16,
    "legend.frameon": False,
})


def main():
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.2))
    for ax, dp in zip(axes, (40, 80)):
        for tag, label, color in CASES:
            eq = OUTDIR / f"{tag}_dp{dp}" / "eq.data"
            if not eq.exists():
                print(f"  {tag}_dp{dp}: missing, skipped")
                continue
            r = char_ratio(str(eq), bond=0.965)
            n = len(r["s"])
            ax.plot(r["s"], r["C"], "-", color=color, lw=3.4,
                    marker="o", ms=10, mfc="white", mew=2.4,
                    markevery=max(n // 9, 1), label=label, zorder=3)
            print(f"  {tag}_dp{dp}: C_app {r['Cinf']:.2f}")
        ax.set_xlabel("contour separation  $s$")
        ax.grid(True, color="0.88", lw=1.0, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(f"DP {dp}")
    axes[0].set_ylabel(r"$C(s) = \langle R^2(s)\rangle / (s\,b^2)$")
    axes[0].legend(loc="lower right", ncol=2, handlelength=1.6)

    # The one caveat that must travel with the right panel: the long chains
    # have not finished relaxing at this run length, so the decline past the
    # peak is equilibration, not lattice structure.
    axes[1].text(0.97, 0.04,
                 "200k NVT steps — DP 80 still equilibrating",
                 transform=axes[1].transAxes, ha="right", va="bottom",
                 fontsize=16, color="0.45", style="italic")

    fig.suptitle("Characteristic ratio of network strands, "
                 "six lattices, no entanglements", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = OUTDIR / "char12_C.png"
    fig.savefig(out, dpi=170, facecolor="white")
    print(f"\n  {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
