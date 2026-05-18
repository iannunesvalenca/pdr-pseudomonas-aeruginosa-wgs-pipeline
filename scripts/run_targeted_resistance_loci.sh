#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"
SEQKIT="${SEQKIT:-$(command -v seqkit || true)}"
BLASTN="${BLASTN:-$(command -v blastn || true)}"

for tool_var in SAMTOOLS SEQKIT BLASTN; do
  if [[ -z "${!tool_var}" ]]; then
    echo "ERROR: ${tool_var} not set and command not found in PATH." >&2
    exit 1
  fi
done

OUTDIR="$ROOT/results/resistance_loci"
mkdir -p "$OUTDIR"

COORDS="$OUTDIR/targeted_loci_coordinates.tsv"
NT_FASTA="$OUTDIR/targeted_loci.nt.fasta"
NT_STATS="$OUTDIR/targeted_loci.nt_stats.tsv"
PAIRWISE="$OUTDIR/targeted_loci.bc42_vs_bc88_blastn.tsv"
AMRFINDER_TABLE="$OUTDIR/targeted_loci_amrfinder_calls.tsv"

cat > "$COORDS" <<'EOF'
sample_id	locus	contig	start	stop	strand	source	note
BC42	oprD	contig_2	17475	18794	+	AMRFinder_all_mutations/Bakta	OprD coding region; Bakta locus OEDHHK_28400
BC88	oprD	contig_1	5444256	5445575	-	AMRFinder_all_mutations/Bakta	OprD coding region; Bakta locus BEHJIH_25255
BC42	mexR	contig_1	4257971	4258306	+	AMRFinder_all_mutations	MexR aligned coding segment; default AMRFinder reports mexR_M112IfsTer8
BC88	mexR	contig_1	2739034	2739369	+	AMRFinder_all_mutations	MexR aligned coding segment; default AMRFinder reports mexR_M112IfsTer8
BC42	mexZ	contig_1	1249894	1250502	+	AMRFinder_all_mutations	MexZ aligned coding region; variant reported only in all_mutations
BC88	mexZ	contig_1	6767511	6768119	+	AMRFinder_all_mutations	MexZ aligned coding region; variant reported only in all_mutations
EOF

assembly_for_sample() {
  case "$1" in
    BC42) printf '%s\n' "$ROOT/results/polished/medaka/BC42/consensus.fasta" ;;
    BC88) printf '%s\n' "$ROOT/results/polished/medaka/BC88/consensus.fasta" ;;
    *) echo "Unknown sample: $1" >&2; return 1 ;;
  esac
}

extract_region() {
  local sample="$1"
  local locus="$2"
  local contig="$3"
  local start="$4"
  local stop="$5"
  local strand="$6"
  local fasta
  fasta="$(assembly_for_sample "$sample")"
  local region="${contig}:${start}-${stop}"
  local header="${sample}|${locus}|${contig}:${start}-${stop}(${strand})"

  if [[ "$strand" == "-" ]]; then
    "$SAMTOOLS" faidx "$fasta" "$region" |
      "$SEQKIT" seq -t DNA -r -p -w 0 |
      awk -v h="$header" 'BEGIN{print ">" h} !/^>/{print}'
  else
    "$SAMTOOLS" faidx "$fasta" "$region" |
      "$SEQKIT" seq -w 0 |
      awk -v h="$header" 'BEGIN{print ">" h} !/^>/{print}'
  fi
}

{
  tail -n +2 "$COORDS" |
    while IFS=$'\t' read -r sample locus contig start stop strand source note; do
      extract_region "$sample" "$locus" "$contig" "$start" "$stop" "$strand"
    done
} > "$NT_FASTA"

"$SEQKIT" fx2tab -n -l -g "$NT_FASTA" > "$NT_STATS"

"$BLASTN" \
  -task blastn \
  -dust no \
  -query "$NT_FASTA" \
  -subject "$NT_FASTA" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen slen qstart qend sstart send evalue bitscore' |
  awk 'BEGIN{FS=OFS="\t"; print "query","subject","pident","alignment_length","mismatch","gapopen","query_length","subject_length","query_coverage_pct","subject_coverage_pct","evalue","bitscore"}
       $1 != $2 {
         split($1, q, "|");
         split($2, s, "|");
         qcov=100*$4/$7;
         scov=100*$4/$8;
         if (q[1] == "BC42" && s[1] == "BC88" && q[2] == s[2] && qcov >= 90 && scov >= 90) {
           print $1,$2,$3,$4,$5,$6,$7,$8,100*$4/$7,100*$4/$8,$13,$14
         }
       }' > "$PAIRWISE"

{
  printf 'sample_id\tsource_file\tsymbol\tcontig\tstart\tstop\tstrand\ttype\tsubtype\tclass\tsubclass\tmethod\tcoverage_pct\tidentity_pct\treference_accession\treference_name\n'
  for sample in BC42 BC88; do
    awk -v sample="$sample" -v source="${sample}_amrfinder.tsv" 'BEGIN{FS=OFS="\t"}
      NR == 1 {next}
      $7 ~ /^mexR_M112IfsTer8$/ || $7 ~ /^nalC_(G71E|S209R)$/ {
        print sample,source,$7,$3,$4,$5,$6,$10,$11,$12,$13,$14,$17,$18,$20,$21
      }' "results/amr/amrfinder/${sample}_amrfinder.tsv"
    awk -v sample="$sample" -v source="${sample}_all_mutations.tsv" 'BEGIN{FS=OFS="\t"}
      NR == 1 {next}
      $7 ~ /^oprD_/ || $7 ~ /^mexR_M112I$/ || $7 ~ /^mexZ_AGLDSLRS194G$/ || $7 ~ /^nalD_A162E$/ {
        print sample,source,$7,$3,$4,$5,$6,$10,$11,$12,$13,$14,$17,$18,$20,$21
      }' "results/amr/amrfinder/${sample}_all_mutations.tsv"
  done
} > "$AMRFINDER_TABLE"

echo "Done. Wrote:"
echo "  $COORDS"
echo "  $NT_FASTA"
echo "  $NT_STATS"
echo "  $PAIRWISE"
echo "  $AMRFINDER_TABLE"
