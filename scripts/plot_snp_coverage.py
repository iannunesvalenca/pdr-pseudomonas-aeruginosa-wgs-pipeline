#!/usr/bin/env python3
"""Read-depth coverage at the 8 assembly-level SNP candidate positions.

Extracts per-base depth from the AP012280-mapped BAMs using samtools depth,
then draws a 2-kb window centred on each candidate position for both isolates.
The confirmed fixed SNP at pos 503340 is highlighted.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

SAMTOOLS = os.environ.get("SAMTOOLS") or shutil.which("samtools")
if SAMTOOLS is None:
    raise SystemExit("ERROR: samtools not found in PATH; set SAMTOOLS=/path/to/samtools")
BAM_DIR = ROOT / "results" / "snp" / "reference_mapping" / "bam"
BC42_BAM = BAM_DIR / "BC42.AP012280.sorted.bam"
BC88_BAM = BAM_DIR / "BC88.AP012280.sorted.bam"
REF_CHROM = "AP012280.1"
OUTDIR = ROOT / "figures" / "parallel_plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

WINDOW = 2000  # bp either side of the SNP position

SNP_POSITIONS = [
    {"pos": 170245,  "label": "170,245\n(mixed)", "status": "artefact"},
    {"pos": 503340,  "label": "503,340\nhslV S90N", "status": "confirmed"},
    {"pos": 1515951, "label": "1,515,951\n(shared)", "status": "artefact"},
    {"pos": 1564187, "label": "1,564,187\n(shared)", "status": "artefact"},
    {"pos": 1956506, "label": "1,956,506\n(mixed)", "status": "artefact"},
    {"pos": 1956619, "label": "1,956,619\n(mixed)", "status": "artefact"},
    {"pos": 1956631, "label": "1,956,631\n(mixed)", "status": "artefact"},
    {"pos": 6536599, "label": "6,536,599\n(shared)", "status": "artefact"},
]

BC42_COLOR = "#2166ac"
BC88_COLOR = "#d73027"
CONFIRMED_BG = "#fff8f0"
ARTEFACT_BG = "#f9f9f9"


def get_depth(bam: Path, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (positions, depths) arrays from samtools depth."""
    region = f"{chrom}:{max(1, start)}-{end}"
    cmd = [SAMTOOLS, "depth", "-r", region, "-d", "0", str(bam)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    positions, depths = [], []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))
    if not positions:
        return np.array([start, end], dtype=int), np.zeros(2, dtype=int)
    return np.array(positions, dtype=int), np.array(depths, dtype=int)


def fill_gaps(pos_arr: np.ndarray, dep_arr: np.ndarray,
              start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
    """Ensure every position from start to end has a depth value (0 if missing)."""
    full_pos = np.arange(start, end + 1, dtype=int)
    full_dep = np.zeros(len(full_pos), dtype=int)
    idx = np.searchsorted(full_pos, pos_arr)
    valid = (idx < len(full_pos)) & (full_pos[idx] == pos_arr)
    full_dep[idx[valid]] = dep_arr[valid]
    return full_pos, full_dep


def main() -> None:
    n = len(SNP_POSITIONS)
    fig, axes = plt.subplots(2, n, figsize=(18, 6.5),
                              gridspec_kw={"hspace": 0.05, "wspace": 0.35})
    fig.patch.set_facecolor("#fafafa")

    for col, snp in enumerate(SNP_POSITIONS):
        pos = snp["pos"]
        start = pos - WINDOW
        end = pos + WINDOW
        status = snp["status"]
        bg = CONFIRMED_BG if status == "confirmed" else ARTEFACT_BG

        for row, (bam, color, sample) in enumerate([
            (BC42_BAM, BC42_COLOR, "BC42"),
            (BC88_BAM, BC88_COLOR, "BC88"),
        ]):
            ax = axes[row][col]
            ax.set_facecolor(bg)

            pos_arr, dep_arr = get_depth(bam, REF_CHROM, start, end)
            pos_full, dep_full = fill_gaps(pos_arr, dep_arr, start, end)
            offset = pos_full - pos  # centre on SNP

            ax.fill_between(offset, dep_full, color=color, alpha=0.30, linewidth=0)
            ax.plot(offset, dep_full, color=color, linewidth=0.8, alpha=0.85)

            # SNP position marker
            ax.axvline(0, color="#333333", linewidth=0.9, linestyle="--", alpha=0.7)

            mean_d = dep_full.mean()
            ax.axhline(mean_d, color=color, linewidth=0.5, linestyle=":",
                       alpha=0.6, zorder=0)

            # Depth annotation
            ax.text(0.97, 0.93, f"μ={mean_d:.0f}×",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=6.5, color=color)

            if status == "confirmed":
                for spine in ax.spines.values():
                    spine.set_edgecolor("#d73027")
                    spine.set_linewidth(1.2)

            ax.set_xlim(-WINDOW, WINDOW)
            ymax = max(dep_full.max() * 1.12, 10)
            ax.set_ylim(0, ymax)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row == 0:
                title_color = "#d73027" if status == "confirmed" else "#555555"
                ax.set_title(snp["label"], fontsize=7, color=title_color,
                             fontweight="bold" if status == "confirmed" else "normal")
                ax.set_xticks([])
            else:
                ax.set_xticks([-WINDOW, 0, WINDOW])
                ax.set_xticklabels([f"-{WINDOW//1000}kb", "SNP", f"+{WINDOW//1000}kb"],
                                   fontsize=6)

            if col == 0:
                ax.set_ylabel(f"{sample}\ndepth (×)", fontsize=7.5,
                              color=color, labelpad=2)
            else:
                ax.set_yticklabels([])

    # Column header: confirmed SNP badge
    confirmed_col = next(i for i, s in enumerate(SNP_POSITIONS)
                         if s["status"] == "confirmed")
    axes[0][confirmed_col].text(
        0.5, 1.35, "✓ Fixed difference", transform=axes[0][confirmed_col].transAxes,
        ha="center", va="bottom", fontsize=7.5, color="#d73027",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", facecolor=CONFIRMED_BG,
                  edgecolor="#d73027", linewidth=0.8),
    )

    legend_handles = [
        mpatches.Patch(facecolor=BC42_COLOR, alpha=0.6, label="BC42 (pre-treatment)"),
        mpatches.Patch(facecolor=BC88_COLOR, alpha=0.6, label="BC88 (post-treatment)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "Read-depth coverage at 8 assembly-level SNP candidate positions (±2 kb window)\n"
        "Mapped to AP012280.1 (NCGM2S1 reference); dashed line = candidate SNP position",
        fontsize=11, fontweight="bold", y=1.02,
    )

    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUTDIR / f"snp_coverage.{ext}",
                    dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())

    print(f"Saved snp_coverage.{{svg,pdf,png}} → {OUTDIR}")


if __name__ == "__main__":
    main()
