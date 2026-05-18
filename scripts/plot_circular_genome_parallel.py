#!/usr/bin/env python3
"""Improved circular genome overview for BC88 — publication-quality version.

Improvements over the original plot_circular_genome_overview.py:
  - gyrA T83I and crpP added to the labelled AMR ring
  - Legend correctly distinguishes AMR (red), regulatory (dark grey), OprD (teal)
  - Label crowding near the dnaA origin resolved by staggering radius
  - Cleaner colour palette with better contrast ratios
  - Tracks slightly separated for visual clarity
  - Output saved to figures/parallel_plots/
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

FIGDIR = ROOT / "figures" / "parallel_plots"
FIGDIR.mkdir(parents=True, exist_ok=True)

BC88_FASTA  = ROOT / "results" / "polished"     / "medaka"     / "BC88" / "consensus.fasta"
BC88_BAKTA  = ROOT / "results" / "annotation"   / "bakta"      / "BC88" / "BC88.tsv"
INTEGRONS   = (
    ROOT / "results" / "mobile_elements" / "integronfinder"
    / "Results_Integron_Finder_BC42_BC88_polished.prefixed"
    / "BC42_BC88_polished.prefixed.integrons"
)
PAF = ROOT / "results" / "visualization" / "synteny" / "BC42_vs_BC88.minimap2.asm5.paf"

COLORS = {
    "cds_fwd":    "#7fa8c9",   # soft steel blue
    "cds_rev":    "#b0b0b0",   # neutral grey
    "amr":        "#c0392b",   # red — resistance genes
    "regulatory": "#4b5563",   # charcoal — mexR / nalC / nalD
    "oprD":       "#0e7c7b",   # teal — OprD
    "gyrA":       "#e05c28",   # burnt orange — QRDR mutation
    "crpP":       "#c97b22",   # amber — acquired fluoroquinolone resistance
    "integron":   "#7b3294",   # purple
    "bc42_plus":  "#9ecae1",   # light blue
    "bc42_minus": "#fdd0a2",   # peach
    "gc_high":    "#1a7340",   # forest green
    "gc_low":     "#8e2c8e",   # dark purple
    "skew_pos":   "#1a56a0",   # dark blue
    "skew_neg":   "#c94b1a",   # dark orange
    "backbone":   "#2c2c2c",
}

# ── All features to label on the AMR ring ──────────────────────────────────
# (display_label, start_1based, end_1based, strand, color)
LABEL_FEATURES = [
    ("Int1: sul1/aadA6",      739807,  742266, "-", COLORS["amr"]),
    ("nalD",                 1733690, 1734328, "+", COLORS["regulatory"]),
    ("nalC",                 1894062, 1894700, "+", COLORS["regulatory"]),
    ("gyrA T83I",            1281900, 1284668, "-", COLORS["gyrA"]),
    ("blaPDC-35 / aph(3')-IIb", 2313804, 2328159, "+", COLORS["amr"]),
    ("mexR",                 2739034, 2741197, "+", COLORS["regulatory"]),
    ("blaOXA-488",           3313511, 3314296, "-", COLORS["amr"]),
    ("crpP",                 4682150, 4682347, "-", COLORS["crpP"]),
    ("Int2: blaIMP-26",      5084416, 5085153, "-", COLORS["amr"]),
    ("oprD",                 5444253, 5445575, "-", COLORS["oprD"]),
]

GENOME_LEN = 7_038_278
ORIGIN0    = 3_246_268   # dnaA start, 0-based


# ── Geometry helpers ───────────────────────────────────────────────────────

def rotate_pos(pos0: int) -> int:
    return (pos0 - ORIGIN0) % GENOME_LEN


def theta(pos0: float) -> float:
    return 2 * math.pi * pos0 / GENOME_LEN


def rotated_segments(start1: int, end1: int) -> list[tuple[int, int]]:
    s0 = start1 - 1
    e0 = end1
    rs = rotate_pos(s0)
    re = rotate_pos(e0)
    if re >= rs:
        return [(rs, re)]
    return [(rs, GENOME_LEN), (0, re)]


def draw_arc(ax, start1: int, end1: int, bottom: float, height: float,
             color: str, alpha: float = 1.0) -> None:
    for rs, re in rotated_segments(start1, end1):
        if re <= rs:
            continue
        mid   = (rs + re) / 2
        width = re - rs
        ax.bar(
            theta(mid), height,
            width=2 * math.pi * width / GENOME_LEN,
            bottom=bottom, align="center",
            color=color, alpha=alpha,
            edgecolor="none", linewidth=0,
        )


def draw_window_bar(ax, start0: int, end0: int, bottom: float,
                    value: float, color: str) -> None:
    mid   = (start0 + end0) / 2
    width = end0 - start0
    ax.bar(
        theta(mid), abs(value),
        width=2 * math.pi * width / GENOME_LEN,
        bottom=bottom if value >= 0 else bottom + value,
        align="center", color=color, edgecolor=color,
        linewidth=0, alpha=0.82,
    )


# ── Data loaders ───────────────────────────────────────────────────────────

def read_tsv(path: Path) -> list[dict[str, str]]:
    header = None
    rows: list[str] = []
    with path.open(newline="") as fh:
        for line in fh:
            if not line.strip():
                continue
            if header is None:
                if line.startswith("#") and not line.startswith("#Sequence Id"):
                    continue
                header = line
            else:
                if line.startswith("#") and not line.startswith("#Sequence Id"):
                    continue
                rows.append(line)
    if header is None:
        return []
    return list(csv.DictReader([header, *rows], delimiter="\t"))


def read_fasta(path: Path) -> str:
    parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line.upper())
    return "".join(parts)


def gc_windows(seq: str, window: int = 20_000) -> list[dict]:
    rotated = seq[ORIGIN0:] + seq[:ORIGIN0]
    rows = []
    for start in range(0, GENOME_LEN, window):
        chunk = rotated[start : min(start + window, GENOME_LEN)]
        if not chunk:
            continue
        gc_n = chunk.count("G") + chunk.count("C")
        gc   = gc_n / len(chunk)
        skew = 0.0 if gc_n == 0 else (chunk.count("G") - chunk.count("C")) / gc_n
        rows.append({"start": start, "end": start + len(chunk), "gc": gc, "skew": skew})
    return rows


def integron_regions() -> list[dict]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in read_tsv(INTEGRONS):
        if row["ID_replicon"] != "BC88|contig_1":
            continue
        grouped[row["ID_integron"]].append((int(row["pos_beg"]), int(row["pos_end"])))
    out = []
    for iid, coords in sorted(grouped.items()):
        s = min(p for pair in coords for p in pair)
        e = max(p for pair in coords for p in pair)
        out.append({"start": s, "end": e})
    return out


def read_paf(path: Path) -> list[dict]:
    rows = []
    if not path.exists():
        return rows
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or int(f[10]) < 10_000:
                continue
            rows.append({"strand": f[4], "start": int(f[7]), "end": int(f[8])})
    return rows


# ── Label placement ────────────────────────────────────────────────────────

def place_labels(ax, features: list[tuple[str, int, int, str]],
                 tick_bottom: float, label_r_base: float,
                 min_sep_rad: float = 0.18) -> None:
    """Draw tick lines and text labels, staggering radius when crowded."""
    # Compute display midpoints and sort by theta
    items = []
    for label, start, end, color in features:
        disp_mid = rotate_pos((start + end) // 2 - 1)
        th = theta(disp_mid)
        items.append({"label": label, "theta": th, "color": color,
                       "start": start, "end": end})
    items.sort(key=lambda x: x["theta"])

    # Assign radii: bump outward if too close to previous label
    radii: list[float] = []
    for i, item in enumerate(items):
        r = label_r_base
        if radii:
            gap = item["theta"] - items[i - 1]["theta"]
            if gap > math.pi:          # wrap-around: measure the short arc
                gap = 2 * math.pi - gap
            if gap < min_sep_rad:
                # Alternate between two radii
                r = label_r_base + 0.055 if (len(radii) % 2 == 0) else label_r_base
        radii.append(r)

    for item, r in zip(items, radii):
        th    = item["theta"]
        color = item["color"]
        label = item["label"]

        # Tick line from track edge to label
        ax.plot([th, th], [tick_bottom, r - 0.01],
                color=color, linewidth=0.75, alpha=0.9)

        # Horizontal alignment: left on right half, right on left half
        x_side = math.sin(th)
        ha     = "left" if x_side >= 0 else "right"
        ax.text(th, r + 0.005, label,
                ha=ha, va="center", fontsize=7.2, color=color,
                fontweight="semibold")


# ── Main ───────────────────────────────────────────────────────────────────

def main() -> None:
    seq      = read_fasta(BC88_FASTA)
    bk_rows  = read_tsv(BC88_BAKTA)

    fig, ax = plt.subplots(figsize=(10, 10),
                            subplot_kw={"projection": "polar"})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2)
    ax.set_ylim(0, 1.20)
    ax.set_axis_off()
    fig.patch.set_facecolor("#fafafa")

    # ── Backbone ring + Mb ticks ──────────────────────────────────────────
    ax.bar(0, 0.010, width=2 * math.pi, bottom=0.970,
           color=COLORS["backbone"], edgecolor="none", alpha=0.85)
    for mb in range(0, int(math.ceil(GENOME_LEN / 1_000_000)) + 1):
        pos = mb * 1_000_000
        if pos > GENOME_LEN:
            continue
        th = theta(pos)
        ax.plot([th, th], [0.960, 0.995], color=COLORS["backbone"],
                linewidth=0.8, alpha=0.7)
        ax.text(th, 1.055, f"{mb}", ha="center", va="center",
                fontsize=7, color="#444444")

    # ── CDS rings ────────────────────────────────────────────────────────
    for row in bk_rows:
        if row.get("Type") != "cds":
            continue
        color  = COLORS["cds_fwd"] if row["Strand"] == "+" else COLORS["cds_rev"]
        bottom = 0.905 if row["Strand"] == "+" else 0.860
        draw_arc(ax, int(row["Start"]), int(row["Stop"]),
                 bottom, 0.030, color, alpha=0.65)

    # ── AMR / key gene markers ────────────────────────────────────────────
    for label, s, e, strand, color in LABEL_FEATURES:
        draw_arc(ax, s, e, 0.795, 0.052, color, alpha=0.95)

    # ── GC content track ─────────────────────────────────────────────────
    windows  = gc_windows(seq)
    mean_gc  = sum(r["gc"] for r in windows) / len(windows)
    gc_max   = max(abs(r["gc"] - mean_gc) for r in windows) or 1
    skew_max = max(abs(r["skew"])          for r in windows) or 1
    # Baseline rings
    ax.bar(0, 0.002, width=2 * math.pi, bottom=0.705,
           color="#888888", edgecolor="none", alpha=0.50)
    ax.bar(0, 0.002, width=2 * math.pi, bottom=0.580,
           color="#888888", edgecolor="none", alpha=0.50)
    for w in windows:
        gc_val = 0.060 * (w["gc"] - mean_gc) / gc_max
        draw_window_bar(ax, w["start"], w["end"], 0.705,
                        gc_val,
                        COLORS["gc_high"] if gc_val >= 0 else COLORS["gc_low"])
        sk_val = 0.060 * w["skew"] / skew_max
        draw_window_bar(ax, w["start"], w["end"], 0.580,
                        sk_val,
                        COLORS["skew_pos"] if sk_val >= 0 else COLORS["skew_neg"])

    # ── BC42 synteny track ────────────────────────────────────────────────
    for row in read_paf(PAF):
        color = COLORS["bc42_plus"] if row["strand"] == "+" else COLORS["bc42_minus"]
        draw_arc(ax, row["start"] + 1, row["end"],
                 0.450, 0.040, color, alpha=0.80)

    # ── IntegronFinder track ──────────────────────────────────────────────
    for region in integron_regions():
        draw_arc(ax, region["start"], region["end"],
                 0.372, 0.048, COLORS["integron"], alpha=0.88)

    # ── Labels ───────────────────────────────────────────────────────────
    label_data = [(lbl, s, e, col)
                  for lbl, s, e, _strand, col in LABEL_FEATURES]
    place_labels(ax, label_data,
                 tick_bottom=0.848, label_r_base=1.070)

    # ── Centre text ───────────────────────────────────────────────────────
    ax.text(0,  0.115, "BC88", ha="center", va="center",
            fontsize=22, fontweight="bold", color="#1a1a1a")
    ax.text(0,  0.045, f"{GENOME_LEN:,} bp", ha="center", va="center",
            fontsize=8.5, color="#444444")
    ax.text(0, -0.020, "origin: dnaA", ha="center", va="center",
            fontsize=7, color="#888888")
    ax.text(0, -0.070, "outer labels: position from dnaA (Mb)",
            ha="center", va="center", fontsize=6.5, color="#aaaaaa")

    # ── Legend ────────────────────────────────────────────────────────────
    legend_handles = [
        Patch(facecolor=COLORS["cds_fwd"],    label="CDS forward strand"),
        Patch(facecolor=COLORS["cds_rev"],    label="CDS reverse strand"),
        Patch(facecolor=COLORS["amr"],        label="AMR gene (acquired)"),
        Patch(facecolor=COLORS["gyrA"],       label="QRDR mutation (gyrA T83I)"),
        Patch(facecolor=COLORS["crpP"],       label="crpP (fluoroquinolone — BLAST confirmed)"),
        Patch(facecolor=COLORS["regulatory"], label="Regulatory (mexR / nalC / nalD)"),
        Patch(facecolor=COLORS["oprD"],       label="OprD (carbapenem porin)"),
        Patch(facecolor=COLORS["integron"],   label="IntegronFinder region"),
        Patch(facecolor=COLORS["bc42_plus"],  label="BC42 alignment — same orientation"),
        Patch(facecolor=COLORS["bc42_minus"], label="BC42 alignment — opposite orientation"),
        Patch(facecolor=COLORS["gc_high"],    label="GC content above mean"),
        Patch(facecolor=COLORS["gc_low"],     label="GC content below mean"),
        Patch(facecolor=COLORS["skew_pos"],   label="GC skew +"),
        Patch(facecolor=COLORS["skew_neg"],   label="GC skew −"),
    ]
    fig.legend(
        legend_handles,
        [h.get_label() for h in legend_handles],
        loc="lower center", ncol=4, frameon=False,
        fontsize=7.5, bbox_to_anchor=(0.5, -0.02),
        handlelength=1.2, handleheight=0.9,
    )

    fig.suptitle("BC88 Circular Genome Overview — PDR P. aeruginosa ST235",
                 fontsize=13, fontweight="bold", y=0.97)
    fig.subplots_adjust(top=0.92, bottom=0.12, left=0.04, right=0.96)

    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIGDIR / f"circular_genome_improved.{ext}",
                    dpi=300, bbox_inches="tight",
                    facecolor=fig.get_facecolor())

    print(f"Saved circular_genome_improved.{{svg,pdf,png}} → {FIGDIR}")


if __name__ == "__main__":
    main()
