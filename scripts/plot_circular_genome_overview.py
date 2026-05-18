#!/usr/bin/env python3
"""Draw a circular BC88 genome overview with annotation and resistance tracks."""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


FIGDIR = ROOT / "figures"
OUTDIR = ROOT / "results" / "visualization" / "circular_genome"
BC88_FASTA = ROOT / "results" / "polished" / "medaka" / "BC88" / "consensus.fasta"
BC88_BAKTA = ROOT / "results" / "annotation" / "bakta" / "BC88" / "BC88.tsv"
BC88_AMRFINDER = ROOT / "results" / "amr" / "amrfinder" / "BC88_amrfinder.tsv"
CRPP_SUMMARY = ROOT / "results" / "blast" / "crpP" / "BC42_BC88_crpP_summary.tsv"
BAKTA_FOCUSED = ROOT / "results" / "annotation" / "bakta_focused" / "BC42_BC88_bakta_key_features.tsv"
INTEGRONS = (
    ROOT
    / "results"
    / "mobile_elements"
    / "integronfinder"
    / "Results_Integron_Finder_BC42_BC88_polished.prefixed"
    / "BC42_BC88_polished.prefixed.integrons"
)
PAF = ROOT / "results" / "visualization" / "synteny" / "BC42_vs_BC88.minimap2.asm5.paf"


COLORS = {
    "cds_forward": "#6f8fb3",
    "cds_reverse": "#a8a8a8",
    "gc_high": "#238b45",
    "gc_low": "#9e3d82",
    "skew_pos": "#225ea8",
    "skew_neg": "#f16913",
    "bc42_plus": "#9ecae1",
    "bc42_minus": "#fdd0a2",
    "integron": "#7b3294",
    "amr": "#c51b29",
    "regulatory": "#4b5563",
    "oprD": "#0f766e",
}

LABEL_BAKTA = {"oprD", "nalD"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    header: str | None = None
    rows: list[str] = []
    with path.open(newline="") as handle:
        for line in handle:
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


def read_fasta(path: Path) -> tuple[str, str]:
    name = ""
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())
    return name, "".join(seq_parts)


def find_dnaa_origin(bakta_rows: list[dict[str, str]]) -> int:
    for row in bakta_rows:
        if row.get("Type") == "cds" and row.get("Gene") == "dnaA":
            return int(row["Start"]) - 1
    for row in bakta_rows:
        if row.get("Type") == "cds" and "replication initiator protein dnaa" in row.get("Product", "").lower():
            return int(row["Start"]) - 1
    return 0


def rotate_pos(pos0: int, origin0: int, genome_len: int) -> int:
    return (pos0 - origin0) % genome_len


def rotated_segments(start1: int, end1: int, origin0: int, genome_len: int) -> list[tuple[int, int]]:
    start0 = start1 - 1
    end0 = end1
    rstart = rotate_pos(start0, origin0, genome_len)
    rend = rotate_pos(end0, origin0, genome_len)
    if rend >= rstart:
        return [(rstart, rend)]
    return [(rstart, genome_len), (0, rend)]


def theta_for_pos(pos0: float, genome_len: int) -> float:
    return 2 * math.pi * pos0 / genome_len


def draw_segment(
    ax,
    start1: int,
    end1: int,
    origin0: int,
    genome_len: int,
    bottom: float,
    height: float,
    color: str,
    alpha: float = 1.0,
    edgecolor: str | None = None,
    linewidth: float = 0.0,
) -> None:
    for rstart, rend in rotated_segments(start1, end1, origin0, genome_len):
        if rend <= rstart:
            continue
        mid = (rstart + rend) / 2
        width = rend - rstart
        ax.bar(
            theta_for_pos(mid, genome_len),
            height,
            width=2 * math.pi * width / genome_len,
            bottom=bottom,
            align="center",
            color=color,
            alpha=alpha,
            edgecolor=edgecolor or color,
            linewidth=linewidth,
        )


def draw_window_bar(ax, start0: int, end0: int, genome_len: int, bottom: float, value: float, color: str) -> None:
    width = end0 - start0
    mid = (start0 + end0) / 2
    ax.bar(
        theta_for_pos(mid, genome_len),
        abs(value),
        width=2 * math.pi * width / genome_len,
        bottom=bottom if value >= 0 else bottom + value,
        align="center",
        color=color,
        edgecolor=color,
        linewidth=0,
        alpha=0.8,
    )


def gc_windows(seq: str, origin0: int, window: int = 20_000) -> list[dict[str, float]]:
    genome_len = len(seq)
    rotated = seq[origin0:] + seq[:origin0]
    rows = []
    for start in range(0, genome_len, window):
        chunk = rotated[start : min(start + window, genome_len)]
        if not chunk:
            continue
        gc = (chunk.count("G") + chunk.count("C")) / len(chunk)
        gc_count = chunk.count("G") + chunk.count("C")
        skew = 0.0 if gc_count == 0 else (chunk.count("G") - chunk.count("C")) / gc_count
        rows.append({"start": start, "end": start + len(chunk), "gc": gc, "skew": skew})
    return rows


def read_paf(path: Path) -> list[dict[str, str]]:
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or int(fields[10]) < 10_000:
                continue
            rows.append(
                {
                    "query": fields[0],
                    "strand": fields[4],
                    "target": fields[5],
                    "start": fields[7],
                    "end": fields[8],
                    "matches": fields[9],
                    "aln_len": fields[10],
                }
            )
    return rows


def integron_regions(origin0: int, genome_len: int) -> list[dict[str, object]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in read_tsv(INTEGRONS):
        if row["ID_replicon"] != "BC88|contig_1":
            continue
        grouped[row["ID_integron"]].append((int(row["pos_beg"]), int(row["pos_end"])))
    regions = []
    for integron_id, coords in sorted(grouped.items()):
        start = min(pos for pair in coords for pos in pair)
        end = max(pos for pair in coords for pos in pair)
        label = "Integron 1" if integron_id.endswith("01") else "Integron 2"
        regions.append({"label": label, "start": start, "end": end})
    return regions


def normalized_symbol(symbol: str) -> str:
    if symbol.startswith("mexR_"):
        return "mexR"
    if symbol.startswith("nalC_"):
        return "nalC"
    if symbol.startswith("gyrA_"):
        return "gyrA"
    return symbol


def label_features() -> list[dict[str, object]]:
    """Return a limited, clustered label set for manuscript readability."""
    by_symbol: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in read_tsv(BC88_AMRFINDER):
        symbol = normalized_symbol(row["Element symbol"])
        by_symbol[symbol].append((int(row["Start"]), int(row["Stop"])))

    by_bakta: dict[str, tuple[int, int]] = {}
    for row in read_tsv(BAKTA_FOCUSED):
        if row["sample_id"] != "BC88" or row["reported_gene"] not in LABEL_BAKTA:
            continue
        by_bakta[row["reported_gene"]] = (int(row["start"]), int(row["stop"]))

    features: list[dict[str, object]] = []

    def add_cluster(label: str, symbols: list[str], color: str) -> None:
        coords: list[tuple[int, int]] = []
        for symbol in symbols:
            coords.extend(by_symbol.get(symbol, []))
        if not coords:
            return
        features.append(
            {
                "label": label,
                "start": min(x for pair in coords for x in pair),
                "end": max(x for pair in coords for x in pair),
                "color": color,
            }
        )

    add_cluster("Int1: sul1/aadA6", ["sul1", "aadA6"], COLORS["amr"])
    add_cluster("gyrA T83I", ["gyrA"], COLORS["amr"])
    add_cluster("blaPDC-35 / aph(3')-IIb", ["blaPDC-35", "aph(3')-IIb"], COLORS["amr"])
    add_cluster("mexR", ["mexR"], COLORS["regulatory"])
    add_cluster("blaOXA-488", ["blaOXA-488"], COLORS["amr"])
    add_cluster("Int2: blaIMP-26", ["blaIMP-26"], COLORS["amr"])
    add_cluster("nalC", ["nalC"], COLORS["regulatory"])

    if CRPP_SUMMARY.exists():
        for row in read_tsv(CRPP_SUMMARY):
            if row.get("sample") != "BC88" or not row.get("call", "").startswith("present"):
                continue
            start = int(row["subject_start"])
            end = int(row["subject_end"])
            features.append({"label": "crpP", "start": min(start, end), "end": max(start, end), "color": COLORS["amr"]})

    for label, color in [("nalD", COLORS["regulatory"]), ("oprD", COLORS["oprD"])]:
        if label in by_bakta:
            start, end = by_bakta[label]
            features.append({"label": label, "start": start, "end": end, "color": color})

    return sorted(features, key=lambda item: int(item["start"]))


def add_label(ax, genome_len: int, origin0: int, start: int, end: int, label: str, color: str) -> None:
    display_mid = rotate_pos((start + end) // 2 - 1, origin0, genome_len)
    theta = theta_for_pos(display_mid, genome_len)
    text_theta = theta
    text_radius = 1.035
    if theta < 0.15 or theta > (2 * math.pi - 0.15):
        text_theta = theta + 0.08 if theta < math.pi else theta - 0.08
        text_radius = 1.075
    # Position label outside the track. Keep text horizontal for readability.
    x_side = math.sin(text_theta)
    ha = "left" if x_side >= 0 else "right"
    ax.plot([theta, theta], [0.93, 1.01], color=color, linewidth=0.7, alpha=0.9)
    if text_theta != theta:
        ax.plot([theta, text_theta], [1.01, text_radius - 0.01], color=color, linewidth=0.55, alpha=0.85)
    ax.text(
        text_theta,
        text_radius,
        label,
        ha=ha,
        va="center",
        fontsize=7.5,
        color=color,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.15},
    )


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    seq_name, seq = read_fasta(BC88_FASTA)
    genome_len = len(seq)
    bakta_rows = read_tsv(BC88_BAKTA)
    origin0 = find_dnaa_origin(bakta_rows)

    fig, ax = plt.subplots(figsize=(9.6, 9.6), subplot_kw={"projection": "polar"})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2)
    ax.set_ylim(0, 1.16)
    ax.set_axis_off()

    # Genome backbone and coordinate ticks.
    ax.bar(0, 0.012, width=2 * math.pi, bottom=0.975, color="#222222", edgecolor="none")
    for mb_tick in range(0, int(math.ceil(genome_len / 1_000_000)) + 1):
        pos = mb_tick * 1_000_000
        if pos > genome_len:
            continue
        theta = theta_for_pos(pos, genome_len)
        ax.plot([theta, theta], [0.965, 1.0], color="#222222", linewidth=0.8)
        ax.text(theta, 1.055, f"{mb_tick}", ha="center", va="center", fontsize=7)

    # CDS rings.
    for row in bakta_rows:
        if row.get("Type") != "cds":
            continue
        color = COLORS["cds_forward"] if row["Strand"] == "+" else COLORS["cds_reverse"]
        bottom = 0.90 if row["Strand"] == "+" else 0.855
        draw_segment(ax, int(row["Start"]), int(row["Stop"]), origin0, genome_len, bottom, 0.032, color, alpha=0.68)

    # GC content and skew tracks.
    windows = gc_windows(seq, origin0)
    mean_gc = sum(row["gc"] for row in windows) / len(windows)
    gc_max = max(abs(row["gc"] - mean_gc) for row in windows) or 1
    skew_max = max(abs(row["skew"]) for row in windows) or 1
    ax.bar(0, 0.002, width=2 * math.pi, bottom=0.705, color="#666666", edgecolor="none", alpha=0.65)
    ax.bar(0, 0.002, width=2 * math.pi, bottom=0.575, color="#666666", edgecolor="none", alpha=0.65)
    for row in windows:
        gc_value = 0.065 * (row["gc"] - mean_gc) / gc_max
        gc_color = COLORS["gc_high"] if gc_value >= 0 else COLORS["gc_low"]
        draw_window_bar(ax, int(row["start"]), int(row["end"]), genome_len, 0.705, gc_value, gc_color)
        skew_value = 0.065 * row["skew"] / skew_max
        skew_color = COLORS["skew_pos"] if skew_value >= 0 else COLORS["skew_neg"]
        draw_window_bar(ax, int(row["start"]), int(row["end"]), genome_len, 0.575, skew_value, skew_color)

    # BC42 alignment and integron tracks.
    for row in read_paf(PAF):
        color = COLORS["bc42_plus"] if row["strand"] == "+" else COLORS["bc42_minus"]
        draw_segment(ax, int(row["start"]) + 1, int(row["end"]), origin0, genome_len, 0.445, 0.042, color, alpha=0.82)

    for region in integron_regions(origin0, genome_len):
        draw_segment(ax, int(region["start"]), int(region["end"]), origin0, genome_len, 0.365, 0.05, COLORS["integron"], alpha=0.9)

    # Resistance/key locus markers and labels.
    feature_table = OUTDIR / "BC88_circular_genome_labeled_features.tsv"
    with feature_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["label", "start", "end", "display_mid_bp"])
        writer.writeheader()
        for feature in label_features():
            start = int(feature["start"])
            end = int(feature["end"])
            draw_segment(ax, start, end, origin0, genome_len, 0.792, 0.055, str(feature["color"]), alpha=0.95)
            add_label(ax, genome_len, origin0, start, end, str(feature["label"]), str(feature["color"]))
            display_mid = rotate_pos((start + end) // 2 - 1, origin0, genome_len)
            writer.writerow({"label": feature["label"], "start": start, "end": end, "display_mid_bp": display_mid})

    ax.text(0, 0.10, "BC88", ha="center", va="center", fontsize=20, fontweight="bold")
    ax.text(0, 0.035, f"{genome_len:,} bp; origin at dnaA", ha="center", va="center", fontsize=8)
    ax.text(0, -0.025, "outer labels: position from dnaA (Mb)", ha="center", va="center", fontsize=7)

    legend_handles = [
        Patch(facecolor=COLORS["cds_forward"], label="CDS forward"),
        Patch(facecolor=COLORS["cds_reverse"], label="CDS reverse"),
        Patch(facecolor=COLORS["amr"], label="AMR/key genes"),
        Patch(facecolor=COLORS["regulatory"], label="Regulatory mutations"),
        Patch(facecolor=COLORS["oprD"], label="OprD porin"),
        Patch(facecolor=COLORS["integron"], label="IntegronFinder region"),
        Patch(facecolor=COLORS["bc42_plus"], label="BC42 same orientation"),
        Patch(facecolor=COLORS["bc42_minus"], label="BC42 opposite orientation"),
        Patch(facecolor=COLORS["gc_high"], label="GC above mean"),
        Patch(facecolor=COLORS["gc_low"], label="GC below mean"),
        Patch(facecolor=COLORS["skew_pos"], label="GC skew +"),
        Patch(facecolor=COLORS["skew_neg"], label="GC skew -"),
    ]
    fig.legend(legend_handles, [h.get_label() for h in legend_handles], loc="lower center", ncol=5, frameon=False, fontsize=8)
    fig.suptitle("BC88 Circular Genome Overview", y=0.98, fontsize=15, fontweight="bold")
    fig.subplots_adjust(top=0.93, bottom=0.10, left=0.04, right=0.96)

    summary = OUTDIR / "BC88_circular_genome_overview_summary.tsv"
    with summary.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"genome\t{seq_name}\n")
        handle.write(f"length_bp\t{genome_len}\n")
        handle.write(f"display_origin\tBC88 dnaA start\n")
        handle.write(f"display_origin_original_bp_0based\t{origin0}\n")
        handle.write(f"mean_gc_window_fraction\t{mean_gc:.6f}\n")
        handle.write(f"gc_window_bp\t20000\n")
        handle.write(f"feature_table\t{feature_table}\n")

    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIGDIR / f"BC88_circular_genome_overview.{ext}", dpi=300, bbox_inches="tight")

    print("Wrote:")
    print(f"  {FIGDIR / 'BC88_circular_genome_overview.svg'}")
    print(f"  {FIGDIR / 'BC88_circular_genome_overview.pdf'}")
    print(f"  {FIGDIR / 'BC88_circular_genome_overview.png'}")
    print(f"  {summary}")
    print(f"  {feature_table}")


if __name__ == "__main__":
    main()
