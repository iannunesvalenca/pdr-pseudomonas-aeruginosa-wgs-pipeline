#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUTDIR="$ROOT/results/annotation/bakta_focused"
mkdir -p "$OUTDIR"

FEATURES="$OUTDIR/BC42_BC88_bakta_key_features.tsv"
COUNTS="$OUTDIR/BC42_BC88_bakta_key_feature_counts.tsv"

{
  printf 'sample_id\tseq_id\ttype\tstart\tstop\tstrand\tlocus_tag\treported_gene\tbakta_gene\tproduct\tdbxrefs\n'
  for sample in BC42 BC88; do
    awk -v sample="$sample" 'BEGIN {
      FS=OFS="\t";
      split("blaIMP-26 blaPDC-35 blaOXA-488 sul1 intI1 oprD mexA mexB oprM mexC mexD oprJ mexE mexF oprN mexX mexY nalC nalD ampR exoU exoT pcrV popB popD spcS vfr toxA toxR nagZ", a, " ");
      for (i in a) keep[a[i]]=1;
    }
    /^#/ { next }
    {
      gene=$7;
      product=$8;
      reported=gene;
      if (gene in keep) {
        print sample,$1,$2,$3,$4,$5,$6,reported,gene,product,$9;
      } else if (product ~ /Adenylate cyclase ExoY/) {
        reported="exoY";
        print sample,$1,$2,$3,$4,$5,$6,reported,gene,product,$9;
      }
    }' "results/annotation/bakta/$sample/$sample.tsv"
  done
} > "$FEATURES"

{
  printf 'sample_id\treported_gene\tcount\n'
  awk 'BEGIN{FS=OFS="\t"} NR>1 {count[$1 OFS $8]++} END{for (k in count) print k,count[k]}' "$FEATURES" |
    sort -k1,1 -k2,2
} > "$COUNTS"

echo "Done. Wrote:"
echo "  $FEATURES"
echo "  $COUNTS"
