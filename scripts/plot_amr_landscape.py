#!/usr/bin/env python3
"""AMR gene landscape — annotated dot chart grouped by drug class.

Shows all resistance genes present in BC42/BC88 (identical between isolates),
grouped by mechanism/drug class. crpP is included from BLAST results.
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUTDIR = ROOT / "figures" / "parallel_plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

AMR_TSV = ROOT / "results" / "amr" / "amrfinder" / "BC88_amrfinder.tsv"

# Manual curation: gene → (display_label, group, color)
# Groups define vertical lanes. crpP added from BLAST result.
CURATED: dict[str, tuple[str, str, str]] = {
    # Beta-lactams
    "blaIMP-26":    ("blaIMP-26\n(MBL)", "Beta-lactam\nresistance", "#c0392b"),
    "blaPDC-35":    ("blaPDC-35\n(AmpC)", "Beta-lactam\nresistance", "#e74c3c"),
    "blaOXA-488":   ("blaOXA-488\n(OXA)", "Beta-lactam\nresistance", "#e57373"),
    # Fluoroquinolones
    "gyrA":         ("gyrA\n(T83I)", "Fluoroquinolone\nresistance", "#e67e22"),
    "crpP":         ("crpP\n(phosphotransferase)", "Fluoroquinolone\nresistance", "#f39c12"),
    # Efflux / regulatory
    "mexR":         ("mexR\n(frameshift)", "Efflux /\nRegulatory", "#8e44ad"),
    "nalC":         ("nalC\n(G71E/S209R)", "Efflux /\nRegulatory", "#9b59b6"),
    "nalD":         ("nalD", "Efflux /\nRegulatory", "#a569bd"),
    "mexA":         ("mexA", "Efflux /\nRegulatory", "#c39bd3"),
    "mexB":         ("mexB", "Efflux /\nRegulatory", "#c39bd3"),
    "mexE":         ("mexE", "Efflux /\nRegulatory", "#c39bd3"),
    "mexX":         ("mexX", "Efflux /\nRegulatory", "#c39bd3"),
    # OprD
    "oprD":         ("oprD\n(intact)", "Porin", "#16a085"),
    # Aminoglycosides
    "aph(3')-IIb":  ("aph(3')-IIb", "Aminoglycoside\nresistance", "#27ae60"),
    "aac(6')-Ib4":  ("aac(6')-Ib4", "Aminoglycoside\nresistance", "#2ecc71"),
    "aacA37":       ("aacA37", "Aminoglycoside\nresistance", "#58d68d"),
    # Sulfonamide / trimethoprim
    "sul1":         ("sul1", "Sulfonamide /\nTrimethoprim", "#2980b9"),
    "dfrA1":        ("dfrA1", "Sulfonamide /\nTrimethoprim", "#3498db"),
    # Chloramphenicol / fosfomycin
    "catB3":        ("catB3", "Other acquired\nresistance", "#7f8c8d"),
    "catB7":        ("catB7", "Other acquired\nresistance", "#95a5a6"),
    "fosA":         ("fosA", "Other acquired\nresistance", "#bdc3c7"),
    "aadA6":        ("aadA6", "Other acquired\nresistance", "#aab7b8"),
    # Mercury (stress)
    "merA":         ("merA", "Metal stress\n(mercury)", "#f0b27a"),
    "merB":         ("merB", "Metal stress\n(mercury)", "#f0b27a"),
    "merD":         ("merD", "Metal stress\n(mercury)", "#f0b27a"),
    "merE":         ("merE", "Metal stress\n(mercury)", "#f0b27a"),
    "merP":         ("merP", "Metal stress\n(mercury)", "#f0b27a"),
    "merR":         ("merR", "Metal stress\n(mercury)", "#f0b27a"),
    "merT":         ("merT", "Metal stress\n(mercury)", "#f0b27a"),
}

GROUP_ORDER = [
    "Beta-lactam\nresistance",
    "Fluoroquinolone\nresistance",
    "Aminoglycoside\nresistance",
    "Sulfonamide /\nTrimethoprim",
    "Efflux /\nRegulatory",
    "Porin",
    "Other acquired\nresistance",
    "Metal stress\n(mercury)",
]

# Mechanism annotation boxes (horizontal bands)
MECHANISM_BANDS = {
    "Beta-lactam\nresistance":      "#fff0f0",
    "Fluoroquinolone\nresistance":  "#fff8f0",
    "Aminoglycoside\nresistance":   "#f0fff4",
    "Sulfonamide /\nTrimethoprim":  "#f0f4ff",
    "Efflux /\nRegulatory":         "#f8f0ff",
    "Porin":                        "#f0fffd",
    "Other acquired\nresistance":   "#f9f9f9",
    "Metal stress\n(mercury)":      "#fffbf0",
}


def load_amr(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def normalize_symbol(sym: str) -> str:
    for prefix in ("mexR_", "nalC_", "gyrA_", "oprD_"):
        if sym.startswith(prefix):
            return sym.split("_")[0]
    return sym


def main() -> None:
    rows = load_amr(AMR_TSV)
    present_symbols = {normalize_symbol(r["Element symbol"]) for r in rows}
    present_symbols.add("crpP")  # confirmed by BLAST

    # Build plot data: group → list of (label, color, identity)
    # identity from AMRFinder for display
    identity_map: dict[str, float] = {}
    for row in rows:
        sym = normalize_symbol(row["Element symbol"])
        try:
            identity_map[sym] = float(row["% Identity to reference"])
        except (ValueError, KeyError):
            pass
    identity_map["crpP"] = 98.48  # from BLAST

    group_genes: dict[str, list[tuple[str, str, float]]] = {g: [] for g in GROUP_ORDER}
    for sym, (label, group, color) in CURATED.items():
        if sym in present_symbols and group in group_genes:
            ident = identity_map.get(sym, 99.0)
            group_genes[group].append((label, color, ident))

    # Layout
    n_groups = len(GROUP_ORDER)
    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    x_positions: list[float] = []
    x_labels: list[str] = []
    x_groups: list[str] = []
    x_colors: list[str] = []
    x_idents: list[float] = []

    gap = 0.6
    x = 0.0
    group_spans: list[tuple[float, float, str]] = []

    for group in GROUP_ORDER:
        genes = group_genes[group]
        if not genes:
            continue
        x_start = x - 0.5
        for label, color, ident in genes:
            x_positions.append(x)
            x_labels.append(label)
            x_groups.append(group)
            x_colors.append(color)
            x_idents.append(ident)
            x += 1.0
        group_spans.append((x_start, x - 0.5, group))
        x += gap

    # Horizontal group bands
    for x_s, x_e, group in group_spans:
        band_color = MECHANISM_BANDS.get(group, "#f5f5f5")
        ax.axvspan(x_s, x_e, color=band_color, alpha=0.6, zorder=0)

    # Lollipop stems
    for xi, yi in zip(x_positions, x_idents):
        ax.plot([xi, xi], [97, yi], color="#cccccc", linewidth=1.0,
                zorder=1, solid_capstyle="round")

    # Dots sized by identity
    sizes = [max((i - 97) * 28 + 60, 20) for i in x_idents]
    sc = ax.scatter(x_positions, x_idents, s=sizes, c=x_colors,
                    zorder=3, linewidths=0.5, edgecolors="white", alpha=0.93)

    # Identity = 100% line
    ax.axhline(100, color="#999999", linewidth=0.6, linestyle=":", zorder=0)

    # Group labels at top
    for x_s, x_e, group in group_spans:
        mid = (x_s + x_e) / 2
        ax.text(mid, 100.5, group, ha="center", va="bottom", fontsize=7.5,
                fontweight="bold", color="#333333")

    # Gene labels (rotated under dots)
    for xi, label in zip(x_positions, x_labels):
        ax.text(xi, 96.8, label, ha="center", va="top", fontsize=6.5,
                rotation=45, color="#333333")

    ax.set_xlim(-0.8, max(x_positions) + 0.8)
    ax.set_ylim(93.5, 101.8)
    ax.set_ylabel("Identity to reference sequence (%)", fontsize=9)
    ax.set_yticks([94, 96, 98, 100])
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # crpP annotation
    crpp_x = x_positions[x_labels.index("crpP\n(phosphotransferase)")]
    ax.annotate("BLAST only\n(AMRFinder miss)",
                xy=(crpp_x, identity_map["crpP"]),
                xytext=(crpp_x + 1.2, identity_map["crpP"] - 0.8),
                fontsize=6.5, color="#e67e22",
                arrowprops=dict(arrowstyle="-", color="#e67e22",
                                lw=0.8, connectionstyle="arc3,rad=0.2"))

    fig.suptitle(
        "Resistance gene landscape — BC42/BC88 (identical profile in both isolates)",
        fontsize=12, fontweight="bold", y=1.01,
    )
    ax.text(0.5, -0.06,
            "All genes confirmed present in both BC42 (pre-treatment) and BC88 (post-treatment). "
            "Dot size reflects identity to reference sequence.",
            ha="center", va="top", transform=ax.transAxes,
            fontsize=7.5, color="#555555", style="italic")

    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUTDIR / f"amr_landscape.{ext}",
                    dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())

    print(f"Saved amr_landscape.{{svg,pdf,png}} → {OUTDIR}")


if __name__ == "__main__":
    main()
