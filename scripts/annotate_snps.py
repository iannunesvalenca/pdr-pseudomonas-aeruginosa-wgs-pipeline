#!/usr/bin/env python3
"""Annotate assembly-level SNPs against the BC88 Bakta annotation."""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNP_TABLE = ROOT / "results/snp/snippy_ctgs/only_diff_snps.tsv"
READ_SUPPORT = ROOT / "results/snp/snippy_ctgs/snippy_vs_read_support_interpretation.tsv"
BC88_FASTA = ROOT / "results/polished/medaka/BC88/consensus.fasta"
REF_FASTA = ROOT / "references/AP012280.1_NCGM2S1.fasta"
REF_GFF = ROOT / "references/AP012280.1_NCGM2S1.gff3"
REF_GFF_FALLBACK = ROOT / "draft_paper_results/snippy/snippy_bc88_ST235/reference/ref.gff"
BC88_BAKTA_TSV = ROOT / "results/annotation/bakta/BC88/BC88.tsv"
BC42_SNIPPY_TAB = ROOT / "results/snp/snippy_ctgs/snippy_BC42_ST235/snippy_BC42_ST235.tab"
BC88_SNIPPY_TAB = ROOT / "results/snp/snippy_ctgs/snippy_BC88_ST235/snippy_BC88_ST235.tab"
OUTDIR = ROOT / "results/snp/snippy_ctgs"
PAF = OUTDIR / "AP012280_vs_BC88.asm5.paf"
OUTPUT = OUTDIR / "snp_annotation_table.tsv"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    return {key: "".join(value) for key, value in seqs.items()}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1].upper()


def translate(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


def run_minimap2() -> None:
    if PAF.exists() and PAF.stat().st_size > 0:
        return
    minimap2 = shutil.which("minimap2")
    if minimap2 is None:
        sys.exit("ERROR: minimap2 not found in PATH")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with PAF.open("w") as out:
        subprocess.run(
            [minimap2, "-x", "asm5", "-c", str(BC88_FASTA), str(REF_FASTA)],
            check=True,
            stdout=out,
        )


def parse_cigar(tags: list[str]) -> list[tuple[int, str]]:
    cg = next((tag[5:] for tag in tags if tag.startswith("cg:Z:")), "")
    if not cg:
        return []
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cg)]


def map_query_pos_to_target(pos1: int, paf_records: list[dict]) -> tuple[str, int, str] | tuple[None, None, str]:
    pos0 = pos1 - 1
    for rec in paf_records:
        if not (rec["qstart"] <= pos0 < rec["qend"]):
            continue
        qpos = rec["qstart"]
        tpos = rec["tstart"] if rec["strand"] == "+" else rec["tend"] - 1
        for length, op in rec["cigar"]:
            consumes_q = op in "MI=X"
            consumes_t = op in "MD=X"
            if consumes_q and qpos <= pos0 < qpos + length:
                if not consumes_t:
                    return None, None, "position falls in query insertion relative to BC88"
                offset = pos0 - qpos
                mapped = tpos + offset if rec["strand"] == "+" else tpos - offset
                return rec["tname"], mapped + 1, "mapped"
            if consumes_q:
                qpos += length
            if consumes_t:
                tpos = tpos + length if rec["strand"] == "+" else tpos - length
        return None, None, "position within alignment but not mappable from CIGAR"
    return None, None, "no AP012280-to-BC88 alignment covers position"


def read_paf(path: Path) -> list[dict]:
    records = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            cigar = parse_cigar(fields[12:])
            if not cigar:
                continue
            records.append(
                {
                    "qname": fields[0],
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tname": fields[5],
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "cigar": cigar,
                }
            )
    return records


def read_bakta_features(path: Path) -> list[dict]:
    features = []
    with path.open() as handle:
        reader = None
        for line in handle:
            if line.startswith("#Sequence Id"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
                break
        if reader is None:
            sys.exit(f"ERROR: no Bakta TSV header found in {path}")
        for row in reader:
            if row["Type"] != "cds":
                continue
            start = int(row["Start"])
            stop = int(row["Stop"])
            if start > stop:
                start, stop = stop, start
            row["_start"] = start
            row["_stop"] = stop
            features.append(row)
    return features


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs = {}
    for item in raw.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def read_reference_gff(path: Path) -> list[dict]:
    features = []
    if not path.exists():
        return features
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gff_attributes(fields[8])
            start = int(fields[3])
            stop = int(fields[4])
            features.append(
                {
                    "seqid": fields[0],
                    "start": min(start, stop),
                    "stop": max(start, stop),
                    "strand": fields[6],
                    "gene": attrs.get("gene") or attrs.get("Name") or "NA",
                    "locus_tag": attrs.get("locus_tag") or attrs.get("ID") or "NA",
                    "product": attrs.get("product") or "NA",
                }
            )
    return features


def same_seqid(a: str, b: str) -> bool:
    return a == b or a.split(".")[0] == b.split(".")[0]


def find_cds(features: list[dict], contig: str, pos1: int) -> dict | None:
    for feature in features:
        if feature["Sequence Id"] == contig and feature["_start"] <= pos1 <= feature["_stop"]:
            return feature
    return None


def find_reference_cds(features: list[dict], chrom: str, pos1: int) -> dict | None:
    for feature in features:
        if same_seqid(feature["seqid"], chrom) and feature["start"] <= pos1 <= feature["stop"]:
            return feature
    return None


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validation_map(path: Path) -> dict[tuple[str, str], dict]:
    return {(row["chrom"], row["pos"]): row for row in read_tsv(path)}


def snippy_type_map(path: Path) -> dict[tuple[str, str], str]:
    if not path.exists():
        return {}
    result = {}
    for row in read_tsv(path):
        result[(row["CHROM"], row["POS"])] = row.get("TYPE", "NA")
    return result


def annotate_coding(feature: dict, genome: dict[str, str], pos1: int, bc42: str, bc88: str) -> dict[str, str]:
    contig = feature["Sequence Id"]
    start = feature["_start"]
    stop = feature["_stop"]
    strand = feature["Strand"]
    genomic = genome[contig][start - 1:stop]
    cds = genomic if strand == "+" else revcomp(genomic)
    offset = pos1 - start if strand == "+" else stop - pos1
    codon_start = (offset // 3) * 3
    codon_index = offset % 3
    bc88_codon = cds[codon_start:codon_start + 3]

    bc42_base = bc42.upper()
    bc88_base = bc88.upper()
    if strand == "-":
        bc42_base = revcomp(bc42_base)
        bc88_base = revcomp(bc88_base)

    bc42_codon = list(bc88_codon)
    if len(bc42_codon) == 3:
        bc42_codon[codon_index] = bc42_base
    bc42_codon = "".join(bc42_codon)

    aa42 = translate(bc42_codon)
    aa88 = translate(bc88_codon)
    codon_number = (offset // 3) + 1
    effect = "synonymous" if aa42 == aa88 else "nonsynonymous"
    if aa42 == "*" and aa88 != "*":
        effect = "stop_lost_in_BC88"
    elif aa42 != "*" and aa88 == "*":
        effect = "stop_gained_in_BC88"

    return {
        "gene": feature.get("Gene", "") or "NA",
        "locus_tag": feature.get("Locus Tag", "") or "NA",
        "product": feature.get("Product", "") or "NA",
        "feature_type": "CDS",
        "strand": strand,
        "codon_number": str(codon_number),
        "bc42_codon": bc42_codon,
        "bc88_codon": bc88_codon,
        "bc42_aa": aa42,
        "bc88_aa": aa88,
        "aa_change_BC42_to_BC88": f"{aa42}{codon_number}{aa88}",
        "effect": effect,
    }


def matching_fasta_key(genome: dict[str, str], seqid: str) -> str | None:
    for key in genome:
        if same_seqid(key, seqid):
            return key
    return None


def annotate_reference_coding(feature: dict, genome: dict[str, str], pos1: int, bc42: str, bc88: str) -> dict[str, str]:
    seqid = matching_fasta_key(genome, feature["seqid"])
    if not seqid:
        return {}
    start = feature["start"]
    stop = feature["stop"]
    strand = feature["strand"]
    genomic = genome[seqid][start - 1:stop]
    cds = genomic if strand == "+" else revcomp(genomic)
    offset = pos1 - start if strand == "+" else stop - pos1
    codon_start = (offset // 3) * 3
    codon_index = offset % 3
    reference_codon = cds[codon_start:codon_start + 3]

    bc42_base = bc42.upper()
    bc88_base = bc88.upper()
    if strand == "-":
        bc42_base = revcomp(bc42_base)
        bc88_base = revcomp(bc88_base)

    bc42_codon = list(reference_codon)
    bc88_codon = list(reference_codon)
    if len(reference_codon) == 3:
        bc42_codon[codon_index] = bc42_base
        bc88_codon[codon_index] = bc88_base
    bc42_codon = "".join(bc42_codon)
    bc88_codon = "".join(bc88_codon)

    aa42 = translate(bc42_codon)
    aa88 = translate(bc88_codon)
    codon_number = (offset // 3) + 1
    effect = "synonymous" if aa42 == aa88 else "nonsynonymous"
    if aa42 == "*" and aa88 != "*":
        effect = "stop_lost_in_BC88"
    elif aa42 != "*" and aa88 == "*":
        effect = "stop_gained_in_BC88"

    return {
        "reference_feature_type": "CDS",
        "reference_gene": feature["gene"],
        "reference_locus_tag": feature["locus_tag"],
        "reference_product": feature["product"],
        "reference_strand": strand,
        "reference_codon_number": str(codon_number),
        "reference_BC42_codon": bc42_codon,
        "reference_BC88_codon": bc88_codon,
        "reference_BC42_aa": aa42,
        "reference_BC88_aa": aa88,
        "reference_aa_change_BC42_to_BC88": f"{aa42}{codon_number}{aa88}",
        "reference_effect": effect,
    }


def main() -> None:
    for required in [SNP_TABLE, READ_SUPPORT, BC88_FASTA, REF_FASTA, BC88_BAKTA_TSV]:
        if not required.exists():
            sys.exit(f"ERROR: required input not found: {required}")

    ref_gff = REF_GFF if REF_GFF.exists() else REF_GFF_FALLBACK
    if not ref_gff.exists():
        sys.exit(f"ERROR: AP012280 reference GFF not found: {REF_GFF} or {REF_GFF_FALLBACK}")

    run_minimap2()
    paf = read_paf(PAF)
    genome = read_fasta(BC88_FASTA)
    reference_genome = read_fasta(REF_FASTA)
    features = read_bakta_features(BC88_BAKTA_TSV)
    reference_features = read_reference_gff(ref_gff)
    validations = validation_map(READ_SUPPORT)
    bc42_snippy_types = snippy_type_map(BC42_SNIPPY_TAB)
    bc88_snippy_types = snippy_type_map(BC88_SNIPPY_TAB)
    snps = read_tsv(SNP_TABLE)

    fieldnames = [
        "chrom",
        "AP012280_pos",
        "ref_allele",
        "BC42_allele",
        "BC88_allele",
        "BC42_snippy_variant_type",
        "BC88_snippy_variant_type",
        "reference_annotation_source",
        "reference_feature_type",
        "reference_gene",
        "reference_locus_tag",
        "reference_product",
        "reference_strand",
        "reference_codon_number",
        "reference_BC42_codon",
        "reference_BC88_codon",
        "reference_BC42_aa",
        "reference_BC88_aa",
        "reference_aa_change_BC42_to_BC88",
        "reference_effect",
        "annotation_caveat",
        "BC88_contig",
        "BC88_pos",
        "BC88_genome_base",
        "mapping_status",
        "feature_type",
        "gene",
        "locus_tag",
        "product",
        "strand",
        "codon_number",
        "BC42_codon",
        "BC88_codon",
        "BC42_aa",
        "BC88_aa",
        "aa_change_BC42_to_BC88",
        "effect",
        "read_validation_status",
        "BC42_depth",
        "BC42_allele_depths",
        "BC88_depth",
        "BC88_allele_depths",
    ]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for snp in snps:
            chrom = snp["CHR"]
            pos = int(snp["POS"])
            bc42 = snp["snippy_BC42_ST235"]
            bc88 = snp["snippy_BC88_ST235"]
            contig, mapped_pos, status = map_query_pos_to_target(pos, paf)
            validation = validations.get((chrom, str(pos)), {})
            reference_feature = find_reference_cds(reference_features, chrom, pos)
            row = {
                "chrom": chrom,
                "AP012280_pos": str(pos),
                "ref_allele": snp["REF"],
                "BC42_allele": bc42,
                "BC88_allele": bc88,
                "BC42_snippy_variant_type": bc42_snippy_types.get((chrom, str(pos)), "none"),
                "BC88_snippy_variant_type": bc88_snippy_types.get((chrom, str(pos)), "none"),
                "reference_annotation_source": str(ref_gff.relative_to(ROOT)),
                "reference_feature_type": "intergenic",
                "reference_gene": "NA",
                "reference_locus_tag": "NA",
                "reference_product": "NA",
                "reference_strand": "NA",
                "reference_codon_number": "NA",
                "reference_BC42_codon": "NA",
                "reference_BC88_codon": "NA",
                "reference_BC42_aa": "NA",
                "reference_BC88_aa": "NA",
                "reference_aa_change_BC42_to_BC88": "NA",
                "reference_effect": "intergenic",
                "annotation_caveat": "none",
                "BC88_contig": contig or "NA",
                "BC88_pos": str(mapped_pos) if mapped_pos else "NA",
                "BC88_genome_base": "NA",
                "mapping_status": status,
                "feature_type": "unmapped",
                "gene": "NA",
                "locus_tag": "NA",
                "product": "NA",
                "strand": "NA",
                "codon_number": "NA",
                "BC42_codon": "NA",
                "BC88_codon": "NA",
                "BC42_aa": "NA",
                "BC88_aa": "NA",
                "aa_change_BC42_to_BC88": "NA",
                "effect": "unmapped",
                "read_validation_status": validation.get("read_backed_interpretation", "NA"),
                "BC42_depth": validation.get("BC42_DP", "NA"),
                "BC42_allele_depths": validation.get("BC42_AD", "NA"),
                "BC88_depth": validation.get("BC88_DP", "NA"),
                "BC88_allele_depths": validation.get("BC88_AD", "NA"),
            }
            if reference_feature:
                row.update(
                    annotate_reference_coding(reference_feature, reference_genome, pos, bc42, bc88)
                )
            variant_types = {row["BC42_snippy_variant_type"], row["BC88_snippy_variant_type"]} - {"none", "snp"}
            if variant_types:
                row["annotation_caveat"] = "reference effect uses the single core-SNP matrix allele; full Snippy variant type: " + ",".join(sorted(variant_types))
            if contig and mapped_pos:
                base = genome[contig][mapped_pos - 1].upper()
                row["BC88_genome_base"] = base
                feature = find_cds(features, contig, mapped_pos)
                if feature:
                    coding = annotate_coding(feature, genome, mapped_pos, bc42, bc88)
                    row.update(
                        {
                            "feature_type": coding["feature_type"],
                            "gene": coding["gene"],
                            "locus_tag": coding["locus_tag"],
                            "product": coding["product"],
                            "strand": coding["strand"],
                            "codon_number": coding["codon_number"],
                            "BC42_codon": coding["bc42_codon"],
                            "BC88_codon": coding["bc88_codon"],
                            "BC42_aa": coding["bc42_aa"],
                            "BC88_aa": coding["bc88_aa"],
                            "aa_change_BC42_to_BC88": coding["aa_change_BC42_to_BC88"],
                            "effect": coding["effect"],
                        }
                    )
                else:
                    row["feature_type"] = "intergenic"
                    row["effect"] = "intergenic"
            writer.writerow(row)

    print(f"SNP annotation table: {OUTPUT}")


if __name__ == "__main__":
    main()
