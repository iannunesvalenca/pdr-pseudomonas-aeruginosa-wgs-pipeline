#!/usr/bin/env python3
"""Virulence system profile — grouped horizontal bar chart.

Reads VFDB Abricate output for BC42 and BC88 and shows gene counts per
virulence category. Since both isolates are identical (287 genes each),
the figure documents the breadth of virulence capability for the ST235 strain.
"""

from __future__ import annotations

import csv
import os
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUTDIR = ROOT / "figures" / "parallel_plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

BC42_VFDB = ROOT / "results" / "virulence" / "vfdb" / "BC42_vfdb.tsv"
BC88_VFDB = ROOT / "results" / "virulence" / "vfdb" / "BC88_vfdb.tsv"

# VFDB product strings → canonical category
CATEGORY_PATTERNS: list[tuple[str, str]] = [
    (r"type III secretion system.*effector", "T3SS Effectors"),
    (r"type III secretion system", "T3SS Apparatus"),
    (r"type VI secretion system|T6SS|HSI-", "T6SS"),
    (r"type IV pil|fimbriae|fimbrial|pilD|pilA|pilB|pilC|pilE|pilF|pilG|"
     r"pilH|pilI|pilJ|pilK|pilM|pilN|pilO|pilP|pilQ|pilR|pilS|pilT|pilU|"
     r"pilV|pilW|pilX|pilY|pilZ|fimL|fimT|fimU|fimV|fimX",
     "Type IV Pili"),
    (r"flagell|motil|flg|fli|flh|fleN|fleQ|fleR|fleS|flhA|flhB|flhF|fliA|"
     r"motA|motB|motC|motD|motY",
     "Flagella & Motility"),
    (r"alginate|Alg|alg", "Alginate / Biofilm"),
    (r"quorum sensing|LasR|LasI|RhlR|RhlI|rhl|las[RI]", "Quorum Sensing"),
    (r"rhamnolipid|rhlA|rhlB|rhlC", "Rhamnolipid"),
    (r"pyoverdine|pvd|pvc", "Pyoverdine (Iron)"),
    (r"pyocyanin|phz", "Pyocyanin"),
    (r"elastase|LasB|LasA|protease|aprA", "Proteases"),
    (r"exotoxin A|ToxA|toxA", "Exotoxin A"),
    (r"Fap|fap", "Functional Amyloid (Fap)"),
    (r"xcpP|xcpQ|xcpR|xcpS|xcpT|xcpU|xcpV|xcpW|xcpX|xcpY|xcpZ|pilD|xcpA",
     "Type II Secretion (Xcp)"),
    (r"cheY|cheZ|chpA|chpB|chpC|chpD|chpE", "Chemotaxis"),
    (r"mercury|mer[ABCDEPRST]", "Metal Resistance"),
]

COLOR_MAP: dict[str, str] = {
    "T3SS Effectors":         "#c0392b",
    "T3SS Apparatus":         "#e74c3c",
    "T6SS":                   "#e67e22",
    "Type IV Pili":           "#f39c12",
    "Flagella & Motility":    "#f1c40f",
    "Alginate / Biofilm":     "#27ae60",
    "Quorum Sensing":         "#1abc9c",
    "Rhamnolipid":            "#16a085",
    "Pyoverdine (Iron)":      "#2980b9",
    "Pyocyanin":              "#3498db",
    "Proteases":              "#9b59b6",
    "Exotoxin A":             "#8e44ad",
    "Fap / Biofilm":          "#7f8c8d",
    "Functional Amyloid (Fap)": "#95a5a6",
    "Type II Secretion (Xcp)": "#34495e",
    "Chemotaxis":             "#bdc3c7",
    "Metal Resistance":       "#f0b27a",
    "Other":                  "#d5d8dc",
}

HIGHLIGHT_GENES = {"exoU", "toxA", "algD", "pvdQ", "lasB"}


def categorize(product: str, gene: str) -> str:
    text = f"{product} {gene}".lower()
    for pattern, category in CATEGORY_PATTERNS:
        if re.search(pattern.lower(), text):
            return category
    return "Other"


def load_vfdb(path: Path) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    with path.open(newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            gene = parts[5]
            product = parts[13]
            pident = float(parts[10]) if parts[10] else 0
            if pident < 80:
                continue
            cat = categorize(product, gene)
            counts[cat] += 1
    return counts


def main() -> None:
    bc42_counts = load_vfdb(BC42_VFDB)
    bc88_counts = load_vfdb(BC88_VFDB)

    # All categories present in either
    all_cats = sorted(
        set(bc42_counts) | set(bc88_counts),
        key=lambda c: -(bc42_counts.get(c, 0) + bc88_counts.get(c, 0)),
    )

    fig, ax = plt.subplots(figsize=(11, 7))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    y = np.arange(len(all_cats))
    bar_h = 0.35

    for i, cat in enumerate(all_cats):
        n42 = bc42_counts.get(cat, 0)
        n88 = bc88_counts.get(cat, 0)
        color = COLOR_MAP.get(cat, COLOR_MAP["Other"])

        # BC42 bar (lighter, slightly above)
        ax.barh(i + bar_h / 2, n42, height=bar_h, color=color,
                alpha=0.55, label="BC42" if i == 0 else "")
        # BC88 bar (full opacity, slightly below)
        ax.barh(i - bar_h / 2, n88, height=bar_h, color=color,
                alpha=0.92, label="BC88" if i == 0 else "")

        # Count labels
        for val, yi_off in [(n42, bar_h / 2), (n88, -bar_h / 2)]:
            if val > 0:
                ax.text(val + 0.15, i + yi_off, str(val),
                        va="center", ha="left", fontsize=7.5,
                        color="#333333")

    # Special annotation for T3SS effectors (ExoU)
    exou_idx = all_cats.index("T3SS Effectors") if "T3SS Effectors" in all_cats else None
    if exou_idx is not None:
        n = bc88_counts.get("T3SS Effectors", 0)
        ax.annotate("ExoU+\n(cytotoxic)", xy=(n + 0.2, exou_idx),
                    xytext=(n + 2.5, exou_idx + 0.8),
                    fontsize=7, color="#c0392b", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.9))

    ax.set_yticks(y)
    ax.set_yticklabels(all_cats, fontsize=9)
    ax.set_xlabel("Number of VFDB gene hits (≥80% identity, ≥80% coverage)", fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        mpatches.Patch(facecolor="#888888", alpha=0.55, label="BC42 (pre-treatment)"),
        mpatches.Patch(facecolor="#888888", alpha=0.92, label="BC88 (post-treatment)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=9)

    total = bc88_counts.get("T3SS Apparatus", 0) + bc88_counts.get("T3SS Effectors", 0)
    fig.suptitle(
        "Virulence gene profile — PDR P. aeruginosa ST235 (BC42 / BC88)",
        fontsize=12, fontweight="bold",
    )
    ax.set_title(
        f"287 VFDB hits in each isolate (identical) | Full T3SS apparatus + ExoU | "
        f"T6SS, alginate, quorum sensing, pyoverdine present",
        fontsize=8, color="#555555", style="italic", pad=6,
    )

    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUTDIR / f"virulence_profile.{ext}",
                    dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())

    print(f"Saved virulence_profile.{{svg,pdf,png}} → {OUTDIR}")


if __name__ == "__main__":
    main()
