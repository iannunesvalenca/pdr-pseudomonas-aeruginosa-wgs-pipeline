#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MICROMAMBA="${PROJECT_ROOT}/.tools/bin/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/micromamba-pdr-pa-wgs}"
export MAMBA_ROOT_PREFIX
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"

SAMPLES=("$@")
if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
  SAMPLES=(BC42 BC88)
fi

cd "${PROJECT_ROOT}"
mkdir -p results/qc/reads logs/qc

run_wgs() {
  "${MICROMAMBA}" run -n pdr-pa-wgs bash -lc "$*"
}

get_sample_field() {
  local sample_id="$1"
  local field_name="$2"
  awk -F '\t' -v sample="${sample_id}" -v field="${field_name}" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == field) {
          col = i
        }
      }
      next
    }
    $1 == sample {
      print $col
      found = 1
    }
    END {
      if (!found) {
        exit 1
      }
    }
  ' config/samples.tsv
}

FASTQS=()
for sample in "${SAMPLES[@]}"; do
  raw_fastq="$(get_sample_field "${sample}" fastq)"
  filtered_fastq="$(get_sample_field "${sample}" filtered_fastq)"
  FASTQS+=("${raw_fastq}" "${filtered_fastq}")
done

seqkit_out="results/qc/reads/$(IFS=_; echo "${SAMPLES[*]}")_seqkit_stats.tsv"
summary_out="results/qc/reads/$(IFS=_; echo "${SAMPLES[*]}")_read_qc_summary.tsv"
nanoq_out="results/qc/reads/$(IFS=_; echo "${SAMPLES[*]}")_nanoq_summary.tsv"

run_wgs "seqkit stats -a -T -j 4 ${FASTQS[*]}" > "${seqkit_out}"

awk -F '\t' '
BEGIN {
  OFS = "\t"
  genome = 7000000
  print "sample_id","read_set","num_reads","yield_bp","yield_mb","estimated_depth_x_7Mb","avg_len","max_len","N50","Q20_pct","Q30_pct","avg_qual","GC_pct"
}
NR > 1 {
  sample = ""
  set = ($1 ~ /filtered/) ? "filtered" : "raw"
  if ($1 ~ /barcode42/) {
    sample = "BC42"
  } else if ($1 ~ /barcode88/) {
    sample = "BC88"
  } else if ($1 ~ /barcode01/) {
    sample = "BC01"
  } else if ($1 ~ /barcode02/) {
    sample = "BC02"
  } else if ($1 ~ /barcode03/) {
    sample = "BC03"
  } else if ($1 ~ /barcode41/) {
    sample = "BC41"
  }
  printf "%s\t%s\t%s\t%s\t%.1f\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", sample,set,$4,$5,$5/1000000,$5/genome,$7,$8,$13,$15,$16,$17,$18
}' "${seqkit_out}" > "${summary_out}"

printf 'sample_id\tread_set\treads\tbases\tn50\tmean_length\tmedian_length\tmean_read_quality\tmedian_read_quality\treads_ge_10kb\treads_ge_30kb\treads_q_ge_15\treads_q_ge_20\n' > "${nanoq_out}"

for sample in "${SAMPLES[@]}"; do
  raw_fastq="$(get_sample_field "${sample}" fastq)"
  filtered_fastq="$(get_sample_field "${sample}" filtered_fastq)"
  for read_set in raw filtered; do
    if [[ "${read_set}" == "raw" ]]; then
      fastq="${raw_fastq}"
    else
      fastq="${filtered_fastq}"
    fi
    json="results/qc/reads/${sample}_${read_set}_nanoq.json"
    run_wgs "nanoq -s -j -i ${fastq}" > "${json}"
    jq -r --arg sample "${sample}" --arg read_set "${read_set}" \
      '[ $sample, $read_set, .reads, .bases, .n50, .mean_length, .median_length, .mean_quality, .median_quality, .length_thresholds["10000"], .length_thresholds["30000"], .quality_thresholds["15"], .quality_thresholds["20"] ] | @tsv' \
      "${json}" >> "${nanoq_out}"
  done
done

echo "Wrote:"
echo "  ${seqkit_out}"
echo "  ${summary_out}"
echo "  ${nanoq_out}"
