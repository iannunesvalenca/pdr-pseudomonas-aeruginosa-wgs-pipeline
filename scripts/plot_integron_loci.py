#!/usr/bin/env python3
"""Draw linear comparison figures for the shared AMR/integron loci."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

FIGDIR = ROOT / "figures"
OUTDIR = ROOT / "results" / "visualization" / "integron_loci"


COLORS = {
    "blaIMP-26": "#c51b29",
    "aadA6": "#2b6cb0",
    "sul1": "#2c7fb8",
    "intI1": "#7b3294",
    "aac(6')-Ib4": "#fdae61",
    "aacA37": "#fdae61",
    "catB3": "#31a354",
    "catB7": "#31a354",
    "dfrA1": "#f6c141",
    "qacEdelta1": "#e07b39",
    "qacG2": "#e07b39",
    "attC": "#111111",
    "other": "#d9d9d9",
}

LABEL_GENES = {
    "blaIMP-26",
    "aadA6",
    "sul1",
    "intI1",
    "aac(6')-Ib4",
    "aacA37",
    "catB3",
    "catB7",
    "dfrA1",
    "qacEdelta1",
    "qacG2",
}


PANELS = [
    {
        "title": "A. sul1/aadA6 class 1 integron region",
        "focus": "sul1/aadA6",
        "regions": ["BC42_sul1_aadA6_context", "BC88_sul1_aadA6_context"],
        "reverse_display": {},
        "identity": "100% nucleotide identity across extracted windows",
    },
    {
        "title": "B. blaIMP-26 integron cassette region",
        "focus": "blaIMP-26/aac/cat/dfr/catB7",
        "regions": ["BC42_blaIMP26_context", "BC88_blaIMP26_context"],
        "reverse_display": {"BC88_blaIMP26_context": True},
        "identity": "100% nucleotide identity across extracted windows; BC88 displayed reverse-complemented",
    },
]


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


def strand_to_sign(strand: str) -> int:
    if strand in {"+", "1", "+1"}:
        return 1
    if strand in {"-", "-1"}:
        return -1
    return 1


def transform(start: int, end: int, region: dict[str, str], reverse: bool) -> tuple[float, float]:
    region_start = int(region["start"])
    region_end = int(region["end"])
    if reverse:
        disp_start = region_end - end + 1
        disp_end = region_end - start + 1
    else:
        disp_start = start - region_start + 1
        disp_end = end - region_start + 1
    return disp_start / 1000, disp_end / 1000


def transform_strand(strand: str, reverse: bool) -> int:
    sign = strand_to_sign(strand)
    return -sign if reverse else sign


def arrow_patch(x0: float, x1: float, y: float, strand: int, color: str, alpha: float = 1.0):
    if x1 < x0:
        x0, x1 = x1, x0
    width = max(x1 - x0, 0.001)
    height = 0.24
    head = min(max(width * 0.35, 0.06), 0.45)
    y0 = y - height / 2
    y1 = y + height / 2
    ym = y
    if strand >= 0:
        pts = [
            (x0, y0),
            (x1 - head, y0),
            (x1, ym),
            (x1 - head, y1),
            (x0, y1),
        ]
    else:
        pts = [
            (x1, y0),
            (x0 + head, y0),
            (x0, ym),
            (x0 + head, y1),
            (x1, y1),
        ]
    return Polygon(pts, closed=True, facecolor=color, edgecolor="#333333", linewidth=0.45, alpha=alpha)


def add_track_label(ax, x: float, y: float, label: str):
    ax.text(x, y, label, ha="right", va="center", fontsize=9, fontweight="bold")


def spread_label_positions(xs: list[float], xmin: float, xmax: float, min_gap: float = 2.25) -> list[float]:
    """Spread crowded gene labels while keeping their order."""
    if not xs:
        return []
    adjusted = [min(max(x, xmin), xmax) for x in xs]
    for i in range(1, len(adjusted)):
        adjusted[i] = max(adjusted[i], adjusted[i - 1] + min_gap)
    overflow = adjusted[-1] - xmax
    if overflow > 0:
        adjusted = [x - overflow for x in adjusted]
        adjusted[0] = max(adjusted[0], xmin)
        for i in range(1, len(adjusted)):
            adjusted[i] = max(adjusted[i], adjusted[i - 1] + min_gap)
    return adjusted


def load_bakta_features(sample: str, region: dict[str, str], reverse: bool) -> list[dict[str, object]]:
    path = ROOT / "results" / "annotation" / "bakta" / sample / f"{sample}.tsv"
    seq_id = region["contig"]
    r_start = int(region["start"])
    r_end = int(region["end"])
    features: list[dict[str, object]] = []

    for row in read_tsv(path):
        if row["#Sequence Id"] != seq_id or row["Type"] != "cds":
            continue
        start = int(row["Start"])
        end = int(row["Stop"])
        if end < r_start or start > r_end:
            continue
        x0, x1 = transform(start, end, region, reverse)
        strand = transform_strand(row["Strand"], reverse)
        gene = row["Gene"] or ""
        product = row["Product"] or ""
        features.append(
            {
                "source": "Bakta",
                "gene": gene,
                "label": gene if gene in LABEL_GENES else "",
                "product": product,
                "x0": x0,
                "x1": x1,
                "strand": strand,
                "color": COLORS.get(gene, COLORS["other"]),
                "z": 1 if gene not in LABEL_GENES else 3,
            }
        )
    return features


def load_amr_features(region_id: str, region: dict[str, str], reverse: bool) -> list[dict[str, object]]:
    amr_path = ROOT / "results" / "amr_context" / "AMR_context_gene_hits.tsv"
    abricate_path = ROOT / "results" / "amr" / "abricate" / "BC42_BC88_ncbi_focused.tsv"
    features: list[dict[str, object]] = []
    seen: set[tuple[str, int, int]] = set()
    for row in read_tsv(amr_path):
        if row["region_id"] != region_id:
            continue
        start = int(row["gene_start"])
        end = int(row["gene_end"])
        gene = row["gene"]
        seen.add((gene, start, end))
        x0, x1 = transform(start, end, region, reverse)
        strand = transform_strand(row["strand"], reverse)
        features.append(
            {
                "source": "AMR",
                "gene": gene,
                "label": gene,
                "product": row["resistance"],
                "x0": x0,
                "x1": x1,
                "strand": strand,
                "color": COLORS.get(gene, "#fb9a99"),
                "z": 5,
            }
        )
    for row in read_tsv(abricate_path):
        if row["sample_id"] != region["sample_id"] or row["contig"] != region["contig"]:
            continue
        start = int(row["start"])
        end = int(row["end"])
        if end < int(region["start"]) or start > int(region["end"]):
            continue
        gene = row["gene"]
        if (gene, start, end) in seen:
            continue
        seen.add((gene, start, end))
        x0, x1 = transform(start, end, region, reverse)
        strand = transform_strand(row["strand"], reverse)
        features.append(
            {
                "source": "Abricate",
                "gene": gene,
                "label": gene,
                "product": row["resistance"],
                "x0": x0,
                "x1": x1,
                "strand": strand,
                "color": COLORS.get(gene, "#fb9a99"),
                "z": 5,
            }
        )

    amrfinder_path = ROOT / "results" / "amr" / "amrfinder" / f"{region['sample_id']}_amrfinder.tsv"
    for row in read_tsv(amrfinder_path):
        if row["Name"] != region["sample_id"] or row["Contig id"] != region["contig"]:
            continue
        gene = row["Element symbol"]
        if gene not in {"qacEdelta1", "qacG2"}:
            continue
        start = int(row["Start"])
        end = int(row["Stop"])
        if end < int(region["start"]) or start > int(region["end"]):
            continue
        if (gene, start, end) in seen:
            continue
        seen.add((gene, start, end))
        x0, x1 = transform(start, end, region, reverse)
        strand = transform_strand(row["Strand"], reverse)
        features.append(
            {
                "source": "AMRFinder",
                "gene": gene,
                "label": gene,
                "product": row["Element name"],
                "x0": x0,
                "x1": x1,
                "strand": strand,
                "color": COLORS[gene],
                "z": 5,
            }
        )
    return features


def load_attc_features(region: dict[str, str], reverse: bool) -> list[dict[str, object]]:
    path = (
        ROOT
        / "results"
        / "mobile_elements"
        / "integronfinder"
        / "Results_Integron_Finder_BC42_BC88_polished.prefixed"
        / "BC42_BC88_polished.prefixed.integrons"
    )
    sample = region["sample_id"]
    replicon = f"{sample}|{region['contig']}"
    r_start = int(region["start"])
    r_end = int(region["end"])
    features: list[dict[str, object]] = []
    for row in read_tsv(path):
        if row["ID_replicon"] != replicon or row["type_elt"] != "attC":
            continue
        start = int(row["pos_beg"])
        end = int(row["pos_end"])
        if end < r_start or start > r_end:
            continue
        x0, x1 = transform(start, end, region, reverse)
        features.append(
            {
                "source": "IntegronFinder",
                "gene": "attC",
                "label": "attC",
                "product": row["annotation"],
                "x0": x0,
                "x1": x1,
                "strand": transform_strand(row["strand"], reverse),
                "color": COLORS["attC"],
                "z": 6,
            }
        )
    return features


def deduplicate_features(features: list[dict[str, object]]) -> list[dict[str, object]]:
    """Keep grey Bakta genes but replace same-locus key genes with curated AMR labels."""
    kept: list[dict[str, object]] = []
    for feat in features:
        if feat["source"] == "Bakta" and feat["gene"] in LABEL_GENES and feat["gene"] != "intI1":
            continue
        else:
            kept.append(feat)
    return kept


def draw_panel(ax, panel: dict[str, object], regions_by_id: dict[str, dict[str, str]], feature_rows: list[dict[str, object]]):
    region_ids = panel["regions"]
    max_len_kb = max((int(regions_by_id[r]["end"]) - int(regions_by_id[r]["start"]) + 1) / 1000 for r in region_ids)
    y_positions = [1.0, 0.0]

    ax.set_title(panel["title"], loc="left", fontsize=11, fontweight="bold", pad=10)
    ax.set_xlim(-7, max_len_kb + 2)
    ax.set_ylim(-0.75, 1.85)
    ax.set_yticks([])
    tick_step = 10 if max_len_kb > 45 else 5
    ax.set_xticks(list(range(0, int(max_len_kb) + tick_step, tick_step)))
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.8)

    for track_idx, (region_id, y) in enumerate(zip(region_ids, y_positions)):
        region = regions_by_id[region_id]
        reverse = bool(panel["reverse_display"].get(region_id, False))
        label_track = track_idx == 0
        length_kb = (int(region["end"]) - int(region["start"]) + 1) / 1000
        track_label = f"{region['sample_id']} {region['contig']}"
        if reverse:
            track_label += " (rev)"
        add_track_label(ax, -0.8, y, track_label)
        ax.add_patch(Rectangle((0, y - 0.035), length_kb, 0.07, color="#999999", alpha=0.45, linewidth=0))

        feats = [f for f in feature_rows if f["region_id"] == region_id]
        labels: list[tuple[float, str]] = []
        for feat in sorted(feats, key=lambda f: (f["z"], f["x0"])):
            if feat["gene"] == "attC":
                mid = (float(feat["x0"]) + float(feat["x1"])) / 2
                ax.plot([mid, mid], [y - 0.24, y + 0.24], color=COLORS["attC"], linewidth=1.1, zorder=8)
                ax.scatter([mid], [y + 0.28], marker="v", s=15, color=COLORS["attC"], zorder=8)
                continue
            alpha = 0.55 if feat["gene"] not in LABEL_GENES else 1.0
            ax.add_patch(
                arrow_patch(
                    float(feat["x0"]),
                    float(feat["x1"]),
                    y,
                    int(feat["strand"]),
                    str(feat["color"]),
                    alpha=alpha,
                )
            )
            if feat["label"] and label_track:
                mid = (float(feat["x0"]) + float(feat["x1"])) / 2
                labels.append((mid, str(feat["label"])))
        sorted_labels = sorted(labels)
        label_xs = spread_label_positions([mid for mid, _ in sorted_labels], 0.4, max_len_kb - 0.4)
        for i, ((mid, label), label_x) in enumerate(zip(sorted_labels, label_xs)):
            level = i % 2
            dy = 0.44 + 0.16 * level if y > 0.5 else -0.44 - 0.16 * level
            va = "bottom" if y > 0.5 else "top"
            ax.plot(
                [mid, label_x],
                [y + (0.18 if y > 0.5 else -0.18), y + dy - (0.03 if y > 0.5 else -0.03)],
                color="#666666",
                linewidth=0.45,
                zorder=9,
            )
            ax.text(label_x, y + dy, label, ha="center", va=va, fontsize=6.8, rotation=38)

    ax.text(max_len_kb, -0.52, str(panel["identity"]), ha="right", va="bottom", fontsize=8, color="#333333")


def write_feature_table(rows: list[dict[str, object]]):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    path = OUTDIR / "integron_locus_figure_features.tsv"
    fieldnames = ["region_id", "sample_id", "contig", "display_reversed", "source", "gene", "x0_kb", "x1_kb", "strand", "label", "product"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def make_feature_rows(regions_by_id: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    all_rows: list[dict[str, object]] = []
    for panel in PANELS:
        for region_id in panel["regions"]:
            region = regions_by_id[region_id]
            reverse = bool(panel["reverse_display"].get(region_id, False))
            features: list[dict[str, object]] = []
            features.extend(load_bakta_features(region["sample_id"], region, reverse))
            features.extend(load_amr_features(region_id, region, reverse))
            features.extend(load_attc_features(region, reverse))
            features = deduplicate_features(features)
            for feat in features:
                feat["region_id"] = region_id
                feat["sample_id"] = region["sample_id"]
                feat["contig"] = region["contig"]
                feat["display_reversed"] = "yes" if reverse else "no"
                feat["x0_kb"] = f"{float(feat['x0']):.3f}"
                feat["x1_kb"] = f"{float(feat['x1']):.3f}"
                all_rows.append(feat)
    return all_rows


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    region_rows = read_tsv(ROOT / "config" / "amr_context_regions.tsv")
    regions_by_id = {row["region_id"]: row for row in region_rows}
    feature_rows = make_feature_rows(regions_by_id)
    table_path = write_feature_table(feature_rows)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(13.2, 8.8))
    for ax, panel in zip(axes, PANELS):
        draw_panel(ax, panel, regions_by_id, feature_rows)

    legend_items = [
        ("blaIMP-26", COLORS["blaIMP-26"]),
        ("aadA6", COLORS["aadA6"]),
        ("sul1", COLORS["sul1"]),
        ("intI1", COLORS["intI1"]),
        ("aac genes", COLORS["aacA37"]),
        ("cat genes", COLORS["catB3"]),
        ("dfrA1", COLORS["dfrA1"]),
        ("qac genes", COLORS["qacEdelta1"]),
        ("other CDS", COLORS["other"]),
        ("attC", COLORS["attC"]),
    ]
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="#333333", linewidth=0.4)
        for _, color in legend_items
    ]
    labels = [label for label, _ in legend_items]
    fig.suptitle("Shared AMR integron loci in BC42 and BC88", fontsize=14, fontweight="bold", y=0.97)
    fig.supxlabel("Position in displayed window (kb)", fontsize=9, y=0.075)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=10, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.15, hspace=0.55)

    for ext in ["svg", "pdf", "png"]:
        out = FIGDIR / f"integron_locus_comparison.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")

    print("Wrote:")
    print(f"  {FIGDIR / 'integron_locus_comparison.svg'}")
    print(f"  {FIGDIR / 'integron_locus_comparison.pdf'}")
    print(f"  {FIGDIR / 'integron_locus_comparison.png'}")
    print(f"  {table_path}")


if __name__ == "__main__":
    main()
