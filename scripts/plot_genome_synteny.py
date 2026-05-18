#!/usr/bin/env python3
"""Draw a compact whole-genome synteny overview for BC42 versus BC88.

BC88 is complete and circular, while BC42 is fragmented. The plot is anchored
to BC88 coordinates. Because the BC88 FASTA origin is arbitrary, the coordinate
system is circularly rotated for display when this avoids splitting a BC42
contig at the plot boundary.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle


BC42 = ROOT / "results" / "polished" / "medaka" / "BC42" / "consensus.fasta"
BC88 = ROOT / "results" / "polished" / "medaka" / "BC88" / "consensus.fasta"
OUTDIR = ROOT / "results" / "visualization" / "synteny"
FIGDIR = ROOT / "figures"
PAF = OUTDIR / "BC42_vs_BC88.minimap2.asm5.paf"


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    seq_len = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
    if name is not None:
        lengths[name] = seq_len
    return lengths


def ensure_paf() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if PAF.exists() and PAF.stat().st_size > 0:
        return
    cmd = [
        "minimap2",
        "-x",
        "asm5",
        "--secondary=no",
        "-t",
        os.environ.get("THREADS", "16"),
        str(BC88),
        str(BC42),
    ]
    with PAF.open("w") as handle:
        subprocess.run(cmd, check=True, stdout=handle)


def read_paf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            aln_len = int(parts[10])
            if aln_len < 10_000:
                continue
            rows.append(
                {
                    "query": parts[0],
                    "qlen": int(parts[1]),
                    "qstart": int(parts[2]),
                    "qend": int(parts[3]),
                    "strand": parts[4],
                    "target": parts[5],
                    "tlen": int(parts[6]),
                    "tstart": int(parts[7]),
                    "tend": int(parts[8]),
                    "matches": int(parts[9]),
                    "aln_len": aln_len,
                    "mapq": int(parts[11]),
                }
            )
    return rows


def mb(value: float) -> float:
    return value / 1_000_000


def draw_scale_bar(ax, x0: float, y: float, length_mb: float = 1.0) -> None:
    ax.plot([x0, x0 + length_mb], [y, y], color="#222222", linewidth=1.2)
    ax.plot([x0, x0], [y - 0.025, y + 0.025], color="#222222", linewidth=1.2)
    ax.plot([x0 + length_mb, x0 + length_mb], [y - 0.025, y + 0.025], color="#222222", linewidth=1.2)
    ax.text(x0 + length_mb / 2, y - 0.07, f"{length_mb:g} Mb", ha="center", va="top", fontsize=8)


def choose_rotation_origin(paf_rows: list[dict[str, object]], target_len: int) -> int:
    """Choose a display origin that joins a query contig split at the circular origin."""
    by_query: dict[str, list[dict[str, object]]] = {}
    for row in paf_rows:
        by_query.setdefault(str(row["query"]), []).append(row)

    for rows in by_query.values():
        if len(rows) < 2:
            continue
        early = [row for row in rows if int(row["tstart"]) < 0.15 * target_len]
        late = [row for row in rows if int(row["tstart"]) > 0.65 * target_len]
        if early and late:
            return min(int(row["tstart"]) for row in late)
    return 0


def rotate_pos(pos: int, origin: int, length: int) -> int:
    return (pos - origin) % length


def add_rotated_coords(row: dict[str, object], origin: int, target_len: int) -> dict[str, object]:
    rotated = dict(row)
    rstart = rotate_pos(int(row["tstart"]), origin, target_len)
    rend = rotate_pos(int(row["tend"]), origin, target_len)
    if rend < rstart:
        raise ValueError("A block crosses the chosen display origin; choose a block boundary as origin.")
    rotated["rstart"] = rstart
    rotated["rend"] = rend
    return rotated


def block_label(row: dict[str, object]) -> str:
    label = str(row["query"])
    if row["strand"] == "-":
        label += " (rev)"
    return label


def display_gaps(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    gaps: list[dict[str, object]] = []
    ordered = sorted(rows, key=lambda r: int(r["rstart"]))
    for left, right in zip(ordered, ordered[1:]):
        start = int(left["rend"])
        end = int(right["rstart"])
        if end <= start:
            continue
        gaps.append(
            {
                "left_block": str(left["query"]),
                "right_block": str(right["query"]),
                "display_start": start,
                "display_end": end,
                "length_bp": end - start,
            }
        )
    return gaps


def add_ice_label(ax, rows: list[dict[str, object]]) -> None:
    for row in rows:
        if row["query"] != "contig_2":
            continue
        start = int(row["rstart"])
        end = int(row["rend"])
        x = mb(start + 0.58 * (end - start))
        ax.plot([x, x], [0.35, 0.53], color="#8a4f20", linewidth=0.7)
        ax.text(
            x,
            0.57,
            "blaIMP-26 ICE\ncontig_2",
            ha="center",
            va="center",
            fontsize=7,
            color="#7a3f13",
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.80, "pad": 0.9},
        )
        return


def merged_label_groups(rows: list[dict[str, object]], max_gap_bp: int = 10_000) -> list[dict[str, object]]:
    groups: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda r: (int(r["rstart"]), int(r["rend"]))):
        label = block_label(row)
        start = int(row["rstart"])
        end = int(row["rend"])
        if groups and groups[-1]["label"] == label and start - int(groups[-1]["end"]) <= max_gap_bp:
            groups[-1]["end"] = max(int(groups[-1]["end"]), end)
            continue
        groups.append({"label": label, "start": start, "end": end, "query": str(row["query"])})
    return groups


def label_for_group(group: dict[str, object]) -> str:
    label = str(group["label"])
    width_bp = int(group["end"]) - int(group["start"])
    if width_bp < 100_000:
        label += f"\n{width_bp / 1000:.0f} kb"
    return label


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ensure_paf()

    bc42_lengths = read_fasta_lengths(BC42)
    bc88_lengths = read_fasta_lengths(BC88)
    paf_rows = read_paf(PAF)

    target_name = next(iter(bc88_lengths))
    target_len = bc88_lengths[target_name]
    rotation_origin = choose_rotation_origin(paf_rows, target_len)
    rotated_rows = [add_rotated_coords(row, rotation_origin, target_len) for row in paf_rows]
    ordered_blocks = sorted(rotated_rows, key=lambda row: (int(row["rstart"]), int(row["rend"])))

    fig, ax = plt.subplots(figsize=(12.2, 4.6))
    y_top = 1.0
    y_bottom = 0.22
    track_h = 0.08

    ax.add_patch(Rectangle((0, y_top - track_h / 2), mb(target_len), track_h, facecolor="#444444", edgecolor="none"))
    top_label = "BC88 complete assembly"
    if rotation_origin:
        top_label += " (rotated)"
    ax.text(-0.05, y_top, top_label, ha="right", va="center", fontsize=9, fontweight="bold")
    ax.text(mb(target_len) + 0.04, y_top, f"{mb(target_len):.2f} Mb", ha="left", va="center", fontsize=8)

    ax.text(-0.05, y_bottom, "BC42 aligned blocks", ha="right", va="center", fontsize=9, fontweight="bold")

    block_table = OUTDIR / "BC42_vs_BC88_synteny_blocks.tsv"
    with block_table.open("w") as handle:
        handle.write(
            "display_order\tbc42_contig\tbc42_start\tbc42_end\tbc42_strand\t"
            "bc88_original_start\tbc88_original_end\tbc88_display_start\tbc88_display_end\t"
            "alignment_length\tidentity_fraction\n"
        )
        for idx, row in enumerate(ordered_blocks, start=1):
            identity = int(row["matches"]) / int(row["aln_len"])
            handle.write(
                f"{idx}\t{row['query']}\t{row['qstart']}\t{row['qend']}\t{row['strand']}\t"
                f"{row['tstart']}\t{row['tend']}\t{row['rstart']}\t{row['rend']}\t"
                f"{row['aln_len']}\t{identity:.6f}\n"
            )

    for idx, row in enumerate(ordered_blocks):
        t0 = mb(int(row["rstart"]))
        t1 = mb(int(row["rend"]))
        color = "#5b8fd6" if row["strand"] == "+" else "#d68f5b"
        poly = Polygon(
            [(t0, y_top - track_h / 2), (t1, y_top - track_h / 2), (t1, y_bottom + track_h / 2), (t0, y_bottom + track_h / 2)],
            closed=True,
            facecolor=color,
            edgecolor="none",
            alpha=0.28,
        )
        ax.add_patch(poly)
        ax.add_patch(Rectangle((t0, y_bottom - track_h / 2), t1 - t0, track_h, facecolor="#777777", edgecolor="none"))

    for idx, group in enumerate(merged_label_groups(ordered_blocks)):
        label = label_for_group(group)
        t0 = mb(int(group["start"]))
        t1 = mb(int(group["end"]))
        block_width = t1 - t0
        mid = (t0 + t1) / 2
        if block_width >= 0.25:
            ax.text(mid, y_bottom - 0.12, label, ha="center", va="top", fontsize=7)
        else:
            y_label = y_bottom + 0.19 + 0.10 * (idx % 2)
            ax.plot([mid, mid], [y_bottom + track_h / 2, y_label - 0.02], color="#555555", linewidth=1.0)
            ax.text(
                mid,
                y_label,
                label,
                ha="center",
                va="bottom",
                fontsize=7.2,
                rotation=32,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.7},
            )

    add_ice_label(ax, ordered_blocks)
    gaps = display_gaps(ordered_blocks)

    max_x = mb(target_len)
    ax.set_xlim(-0.25, max_x + 0.25)
    ax.set_ylim(-0.05, 1.28)
    ax.set_yticks([])
    xlabel = "Rotated BC88 genome position (Mb)" if rotation_origin else "BC88 genome position (Mb)"
    ax.set_xlabel(xlabel, fontsize=9)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.8)
    ax.set_title("Whole-genome synteny: BC42 alignment blocks on BC88", loc="left", fontsize=12, fontweight="bold")
    identities = [int(row["matches"]) / int(row["aln_len"]) for row in paf_rows]
    max_gap = max((int(gap["length_bp"]) for gap in gaps), default=0)
    ax.text(
        0,
        1.17,
        f"Alignment block identity: {min(identities) * 100:.2f}-{max(identities) * 100:.2f}%",
        ha="left",
        va="center",
        fontsize=8,
        color="#333333",
    )
    ax.text(
        0,
        1.09,
        f"White gaps are BC88 intervals without a displayed BC42 block; largest gap {max_gap:,} bp",
        ha="left",
        va="center",
        fontsize=7.2,
        color="#555555",
    )

    ax.add_patch(Rectangle((max_x - 1.95, 1.18), 0.13, 0.06, facecolor="#5b8fd6", alpha=0.45, edgecolor="none"))
    ax.text(max_x - 1.76, 1.21, "same orientation", ha="left", va="center", fontsize=8)
    ax.add_patch(Rectangle((max_x - 0.95, 1.18), 0.13, 0.06, facecolor="#d68f5b", alpha=0.45, edgecolor="none"))
    ax.text(max_x - 0.76, 1.21, "opposite orientation", ha="left", va="center", fontsize=8)
    draw_scale_bar(ax, 0, 0.02)

    summary = OUTDIR / "BC42_vs_BC88_synteny_summary.tsv"
    gap_table = OUTDIR / "BC42_vs_BC88_synteny_display_gaps.tsv"
    with gap_table.open("w") as handle:
        handle.write("left_block\tright_block\tdisplay_start\tdisplay_end\tlength_bp\n")
        for gap in gaps:
            handle.write(
                f"{gap['left_block']}\t{gap['right_block']}\t{gap['display_start']}\t"
                f"{gap['display_end']}\t{gap['length_bp']}\n"
            )
    aligned_bp = sum(int(row["aln_len"]) for row in paf_rows)
    with summary.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write("display_mode\tBC42 alignment blocks ordered by rotated BC88 coordinate\n")
        handle.write(f"bc88_display_origin_original_bp\t{rotation_origin}\n")
        handle.write(f"BC88_length_bp\t{target_len}\n")
        handle.write(f"BC42_total_length_bp\t{sum(bc42_lengths.values())}\n")
        handle.write(f"BC42_contigs\t{len(bc42_lengths)}\n")
        handle.write(f"alignment_blocks_ge_10kb\t{len(paf_rows)}\n")
        handle.write(f"sum_alignment_block_lengths_bp\t{aligned_bp}\n")
        handle.write(f"paf\t{PAF}\n")
        handle.write(f"block_table\t{block_table}\n")
        handle.write(f"gap_table\t{gap_table}\n")

    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIGDIR / f"BC42_BC88_whole_genome_synteny.{ext}", dpi=300, bbox_inches="tight")

    print("Wrote:")
    print(f"  {FIGDIR / 'BC42_BC88_whole_genome_synteny.svg'}")
    print(f"  {FIGDIR / 'BC42_BC88_whole_genome_synteny.pdf'}")
    print(f"  {FIGDIR / 'BC42_BC88_whole_genome_synteny.png'}")
    print(f"  {summary}")


if __name__ == "__main__":
    main()
