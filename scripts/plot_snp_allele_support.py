#!/usr/bin/env python3
"""SNP candidate allele-support panel — redesigned for clarity.

Three visually distinct categories:
  fixed  — BC42 and BC88 carry OPPOSITE alleles (real SNP)
  shared — both isolates carry the SAME non-reference allele (divergence from
            the AP012280 reference, not a BC42 vs BC88 difference)
  mixed  — ~50/50 split in both samples (repeat / alignment noise)
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

OUTDIR = ROOT / "figures" / "parallel_plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

SNP_DATA = [
    {
        "pos": "170,245",
        "label": "putative\ndeoxyribosyl-\ntransferase",
        "bc42": (115, 53), "bc88": (70, 42),
        "ref_allele": "C", "alt_allele": "T",
        "status": "mixed",
    },
    {
        "pos": "503,340",
        "label": "hslV\n(S90N)",
        "bc42": (142, 0), "bc88": (0, 133),
        "ref_allele": "C", "alt_allele": "T",
        "status": "fixed",
    },
    {
        "pos": "1,515,951",
        "label": "hypothetical\nprotein",
        "bc42": (2, 119), "bc88": (3, 73),
        "ref_allele": "G", "alt_allele": "C",
        "status": "shared",
    },
    {
        "pos": "1,564,187",
        "label": "hypothetical\nprotein",
        "bc42": (7, 131), "bc88": (2, 92),
        "ref_allele": "C", "alt_allele": "G",
        "status": "shared",
    },
    {
        "pos": "1,956,506",
        "label": "intergenic",
        "bc42": (86, 65), "bc88": (79, 69),
        "ref_allele": "G", "alt_allele": "C",
        "status": "mixed",
    },
    {
        "pos": "1,956,619",
        "label": "intergenic",
        "bc42": (87, 73), "bc88": (75, 76),
        "ref_allele": "T", "alt_allele": "C",
        "status": "mixed",
    },
    {
        "pos": "1,956,631",
        "label": "intergenic",
        "bc42": (85, 74), "bc88": (73, 77),
        "ref_allele": "G", "alt_allele": "T",
        "status": "mixed",
    },
    {
        "pos": "6,536,599",
        "label": "intergenic",
        "bc42": (21, 117), "bc88": (19, 92),
        "ref_allele": "G", "alt_allele": "A",
        "status": "shared",
    },
]

# Per-status colour palette
STYLE = {
    "fixed": {
        "bc42_ref": "#2166ac",  "bc42_alt": "#acd3f0",   # blue family
        "bc88_ref": "#f4b8b0",  "bc88_alt": "#c0392b",   # red family
        "bg": "#fff7f6",
        "badge_fc": "#c0392b",  "badge_ec": "#922b21",
        "badge_text": "FIXED DIFFERENCE",
        "badge_icon": "✓",
        "spine_color": "#c0392b",
    },
    "shared": {
        "bc42_ref": "#d5e8f7",  "bc42_alt": "#e67e22",   # orange for the shared allele
        "bc88_ref": "#d5e8f7",  "bc88_alt": "#e67e22",
        "bg": "#fffbf5",
        "badge_fc": "#e67e22",  "badge_ec": "#b8620d",
        "badge_text": "SHARED VARIANT",
        "badge_icon": "=",
        "spine_color": "#e67e22",
    },
    "mixed": {
        "bc42_ref": "#aaaaaa",  "bc42_alt": "#dddddd",
        "bc88_ref": "#aaaaaa",  "bc88_alt": "#dddddd",
        "bg": "#f9f9f9",
        "badge_fc": "#888888",  "badge_ec": "#555555",
        "badge_text": "MIXED / NOISE",
        "badge_icon": "~",
        "spine_color": "#aaaaaa",
    },
}


def dominant_allele(ref: int, alt: int, ref_base: str, alt_base: str) -> str:
    return alt_base if alt >= ref else ref_base


def dominant_pct(ref: int, alt: int) -> float:
    total = ref + alt
    return max(ref, alt) / total if total else 0.0


def draw_bar(ax, x, ref, alt, style_ref, style_alt, alpha, bar_w=0.30):
    total = ref + alt
    if total == 0:
        return
    ax.bar(x, ref / total, width=bar_w, bottom=0,
           color=style_ref, alpha=alpha, linewidth=0)
    ax.bar(x, alt / total, width=bar_w, bottom=ref / total,
           color=style_alt, alpha=alpha, linewidth=0)


fig, axes = plt.subplots(1, 8, figsize=(17, 6.0),
                          gridspec_kw={"wspace": 0.60})
fig.patch.set_facecolor("#f7f7f7")

BAR_H   = 1.0
X42     = -0.22
X88     =  0.22

for idx, (snp, ax) in enumerate(zip(SNP_DATA, axes)):
    st     = snp["status"]
    sty    = STYLE[st]
    alpha  = 1.0 if st == "fixed" else 0.82
    ref_a  = snp["ref_allele"]
    alt_a  = snp["alt_allele"]

    ax.set_facecolor(sty["bg"])

    # ── Bars ──────────────────────────────────────────────────────────────
    draw_bar(ax, X42, snp["bc42"][0], snp["bc42"][1],
             sty["bc42_ref"], sty["bc42_alt"], alpha)
    draw_bar(ax, X88, snp["bc88"][0], snp["bc88"][1],
             sty["bc88_ref"], sty["bc88_alt"], alpha)

    # ── Dominant allele letter inside bar ─────────────────────────────────
    for x, (ref, alt) in [(X42, snp["bc42"]), (X88, snp["bc88"])]:
        dom = dominant_allele(ref, alt, ref_a, alt_a)
        pct = dominant_pct(ref, alt)
        # Place letter at bar centre
        ax.text(x, 0.50, dom, ha="center", va="center",
                fontsize=11, fontweight="bold",
                color="white" if pct > 0.65 else "#333333",
                alpha=alpha)
        # Purity % above bar
        ax.text(x, BAR_H + 0.06, f"{pct:.0%}",
                ha="center", va="bottom", fontsize=7.5,
                color="#333333", alpha=alpha)

    # ── For SHARED: bracket + "same allele" annotation ────────────────────
    if st == "shared":
        shared_base = dominant_allele(*snp["bc42"], ref_a, alt_a)
        ax.annotate(
            "", xy=(X88 - 0.01, 1.09), xytext=(X42 + 0.01, 1.09),
            arrowprops=dict(arrowstyle="<->", color="#e67e22",
                            lw=1.2, connectionstyle="arc3,rad=0")
        )
        ax.text(0, 1.14, f"both = {shared_base}",
                ha="center", va="bottom", fontsize=6.8,
                color="#b8620d", fontweight="bold")

    # ── Reference allele footnote ──────────────────────────────────────────
    ax.text(0, -0.12, f"ref: {ref_a}",
            ha="center", va="top", fontsize=6.5, color="#888888")

    # ── Spine styling ─────────────────────────────────────────────────────
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_edgecolor(sty["spine_color"])
        ax.spines[spine].set_linewidth(1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Verdict badge (top of panel) ──────────────────────────────────────
    badge_y = 1.35
    ax.text(0, badge_y,
            f"{sty['badge_icon']}  {sty['badge_text']}",
            ha="center", va="center", fontsize=6.0,
            fontweight="bold", color="white",
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor=sty["badge_fc"],
                edgecolor=sty["badge_ec"],
                linewidth=0.8,
            ),
            transform=ax.transData,
    )

    # ── Axes ──────────────────────────────────────────────────────────────
    ax.set_xlim(-0.58, 0.58)
    ax.set_ylim(-0.30, 1.55)
    ax.set_xticks([X42, X88])
    ax.set_xticklabels(["BC42", "BC88"], fontsize=8)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"] if idx == 0 else [], fontsize=7)
    if idx == 0:
        ax.set_ylabel("Read fraction (ref ↓ / alt ↑)", fontsize=7.5)

    ax.axhline(0.5, color="#cccccc", linewidth=0.6, linestyle="--", zorder=0)

    # Position + gene label as title
    ax.set_title(f"{snp['pos']}\n{snp['label']}",
                 fontsize=6.6, pad=3, color="#222222")

# ── Global legend ──────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(facecolor="#2166ac", label="BC42 — reference allele"),
    mpatches.Patch(facecolor="#c0392b", label="BC88 — alternate allele"),
    mpatches.Patch(facecolor="#e67e22", label="Shared allele (both isolates identical)"),
    mpatches.Patch(facecolor="#aaaaaa", label="Mixed reads (~50/50, ambiguous)"),
    mpatches.Patch(facecolor="#d5e8f7", label="Reference allele (minority)"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=5,
           frameon=False, fontsize=7.8, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(
    "Read-level allele support at 8 assembly-level SNP candidates  ·  "
    "BC42 (pre-treatment) vs BC88 (post-treatment)\n"
    "Letters show dominant allele; % shows read purity for that allele",
    fontsize=10.5, fontweight="bold", y=1.01,
)

for ext in ["svg", "pdf", "png"]:
    fig.savefig(OUTDIR / f"snp_allele_support.{ext}",
                dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())

print(f"Saved snp_allele_support.{{svg,pdf,png}} → {OUTDIR}")
