#!/usr/bin/env bash
# Demo data for the VESPA workflow: two 100 kb windows of HG008 chr20, each holding
# one somatic SV of the GIAB V0.5 draft benchmark, as two contigs.
#   demo_del  chr20:51,746,001-51,846,000  DEL 1368 bp at 51,796,494 (placed path)
#   demo_inv  chr20:38,379,001-38,479,000  INV 1359 bp at 38,428,242 (clipped path)
# Reads: primary alignments contained in a window, from the NIST Revio jan2024 BAMs.
# Phased VCF: the normal.phased.vcf.gz of the workflow (DeepVariant + whatshap phase) run on
# chr20 of the same BAMs. Needs samtools >= 1.14 (endpos in filter expressions).
# Usage: make_demo.sh TUMOR.bam NORMAL.bam chr20.fa normal.phased.vcf.gz OUTDIR
set -euo pipefail
T=$1 N=$2 REF=$3 VCF=$4 O=$5
mkdir -p "$O"
declare -A REG=([demo_del]=51746001 [demo_inv]=38379001)
LEN=100000
: > "$O/ref.fa"
for s in tumor normal; do : > "$O/$s.fq"; done
vcf_body=$(mktemp)
for c in demo_del demo_inv; do
  b=${REG[$c]} e=$(( ${REG[$c]} + LEN - 1 )) r="chr20:${REG[$c]}-$(( ${REG[$c]} + LEN - 1 ))"
  samtools faidx "$REF" "$r" | sed "1s/.*/>$c/" >> "$O/ref.fa"
  for s in tumor normal; do
    bam=$T; [ $s = normal ] && bam=$N
    samtools view -u -F 0x904 -e "pos >= $b && endpos <= $e" "$bam" "$r" \
      | samtools fastq - 2>/dev/null >> "$O/$s.fq"
  done
  bcftools view -H -r "$r" "$VCF" | awk -v OFS='\t' -v c=$c -v o=$(( b - 1 )) '{ $1 = c; $2 -= o; print }' >> "$vcf_body"
done
samtools faidx "$O/ref.fa"
for s in tumor normal; do gzip -9f "$O/$s.fq"; done
{
  echo '##fileformat=VCFv4.2'
  bcftools view -h "$VCF" | grep -E '^##(FORMAT|INFO|FILTER)='
  for c in demo_del demo_inv; do echo "##contig=<ID=$c,length=$LEN>"; done
  bcftools view -h "$VCF" | tail -n 1
  cat "$vcf_body"
} | bgzip > "$O/normal.phased.vcf.gz"
tabix -f -p vcf "$O/normal.phased.vcf.gz"
rm -f "$vcf_body"
