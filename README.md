# VESPA: Variant Extraction from somatic-SPecific strings and draft Assemblies

#### Note: VESPA needs a tumor and its matched normal, both sequenced with PacBio HiFi. It has not been tested on other technologies (e.g., ONT).

---

VESPA finds the sequence of the tumor reads that is absent from the reference and from the germline genome of the patient, and calls somatic SVs from it.
It is based on the notion of sample-specific strings (SFS, or simply _specific strings_): the shortest substrings of a read that do not occur in an index.
Here the index holds the reference and the germline sequence of the patient: a draft assembly of the normal reads (unitigs) and the normal reads the assembly misses.
The SFS of the tumor reads against this index are then specific to the tumor, and they mark candidate somatic SVs.

The steps are:
1. error correction of the tumor and normal reads with a version of hifiasm that preserves somatic variants ([hifiasm v1.0.0-somatic](https://github.com/ldenti/hifiasm/tree/v1.0.0-somatic));
2. index of the reference, the normal unitigs and the normal reads they miss; SFS of the tumor reads against it;
3. calling from the SFS placed on the tumor alignments (clustering and partial-order assembly) and from the SFS at clipped alignments (breakpoint junctions);
4. germline filtering with the normal reads.

VESPA extends [SVDSS](https://github.com/Parsoa/SVDSS) (Denti, Khorsand et al., [Nature Methods 2023](https://doi.org/10.1038/s41592-022-01674-1)), which calls germline SVs from SFS.

## Download and Installation

You can get VESPA in two ways:
* [downloading a static binary](#static-binary)
* [compiling it](#compilation-from-source)

#### Static Binary
A static binary for x86_64 Linux, `vespa_linux_x86-64`, is attached to each [release](https://github.com/dcmonti/VESPA/releases/latest).
It has no runtime dependency.

#### Compilation from Source
To compile VESPA, you need:
* a C++14-compliant compiler (GCC 8.2 or newer)
* make, automake, autoconf
* cmake (>=3.14)
* git
* some other development libraries: zlib, bz2, lzma

To install these dependencies:
```bash
# On a deb-based system:
sudo apt install build-essential autoconf cmake git zlib1g-dev libbz2-dev liblzma-dev
# On a rpm-based system:
sudo dnf install gcc gcc-c++ make automake autoconf cmake git libstdc++-static zlib-devel bzip2-devel xz-devel
```

The following libraries are downloaded and compiled while compiling VESPA, each at a fixed commit:
* [htslib](https://github.com/samtools/htslib) built with [libdeflate](https://github.com/ebiggers/libdeflate) for BAM processing.
* [ksw2](https://github.com/lh3/ksw2) for alignment.
* [ropebwt3](https://github.com/lh3/ropebwt3) for FMD index creation and querying.
* [abPOA](https://github.com/yangao07/abPOA) for POA computation.
* [rapidfuzz](https://github.com/maxbachmann/rapidfuzz-cpp) for string similarity computation.
* [interval-tree](https://github.com/5cript/interval-tree) for variant overlap detection and clustering.
* [spdlog](https://github.com/gabime/spdlog) for logging.

To download and compile VESPA:
```bash
git clone https://github.com/dcmonti/VESPA.git
cd VESPA
mkdir build ; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
This creates the `vespa` binary in the root of the repo.
`static/build-static.sh` builds the static binary instead.

## Usage Guide

The input is the HiFi reads of a tumor and of its matched normal (FASTQ, FASTA or unaligned BAM converted to FASTQ) and a reference genome.
Besides `vespa`, the steps below use [hifiasm v1.0.0-somatic](https://github.com/ldenti/hifiasm/tree/v1.0.0-somatic), [minimap2](https://github.com/lh3/minimap2), [samtools](https://github.com/samtools/samtools), [bcftools](https://github.com/samtools/bcftools), Python 3 with [pysam](https://github.com/pysam-developers/pysam) for the helpers in `scripts/`, and, for haplotagging, [DeepVariant](https://github.com/google/deepvariant) and [whatshap](https://github.com/whatshap/whatshap).

In the guide below the reference is `ref.fa` (the chromosomes to call on) and the reads are `tumor.fq.gz` and `normal.fq.gz`.

### 1. Error correction

```bash
hifiasm -t16 --write-ec --het-cov <N> --bin-only -o tumor tumor.fq.gz    # tumor.ec.fa
hifiasm -t16 --write-ec --het-cov <N> -o normal normal.fq.gz             # normal.ec.fa, normal unitigs
awk '/^S/ {print ">" $2; print $3}' normal.bp.r_utg.gfa > normal.r_utg.fa
```
`--het-cov` keeps the somatic variants supported by few reads from being corrected; for the tumor it can instead be computed from a negative-binomial model (`--tumor-purity`, `--clonality`, `--dispersion`, `--somatic-fnr`).
The options are described in the [README of the fork](https://github.com/ldenti/hifiasm/tree/het_cov#somatic-error-correction).

### 2. Alignment and haplotagging

```bash
for s in tumor normal; do
  minimap2 -ax map-hifi -t16 ref.fa $s.ec.fa | samtools sort -o $s.ec.bam - && samtools index $s.ec.bam
done
```
Haplotagging is optional; with it, VESPA keeps the SFS of the two haplotypes apart and the germline filters use the haplotype of each normal read.
Small variants of the normal are called with DeepVariant on the raw normal alignments and phased with `whatshap phase`; `whatshap haplotag --output-haplotag-list` assigns the raw reads of each sample to a haplotype, and the corrected reads, which keep the name of their raw read, are tagged from that list:
```bash
python scripts/tag_bam_from_list.py tumor.ec.bam tumor.haplotags.tsv tumor.ec.hp.bam --with-ps && samtools index tumor.ec.hp.bam
```

### 3. Index and SFS

```bash
# base index: reference + normal unitigs, homopolymers capped at 2
vespa index -d -m 4G -t16 -p 0 --hpc --hpc-cap 2 ref.fa normal.r_utg.fa -o base.fmd
# normal reads that still produce SFS: germline sequence the unitigs miss
vespa search --index base.fmd --fastx normal.ec.fa --threads 16 --hpc --hpc-cap 2 > normal.sfs
python scripts/extract_missing_reads.py normal.sfs normal.ec.fa missing.fa.gz
# full index: base index + those reads
vespa index -d -m 4G -t16 -p 0 --hpc --hpc-cap 2 -i base.fmd missing.fa.gz -o full.fmd
# SFS of the tumor reads, deduplicated, without SFS longer than 500 bp
vespa search --index full.fmd --fastx tumor.ec.fa --threads 16 --hpc --hpc-cap 2 > tumor.sfs
python scripts/dedup_sfs.py tumor.sfs tumor.dedup.sfs --max-length 500
```
The `index` subcommand is a wrapper around `ropebwt3 build` (it exposes the same interface), with `--hpc` and `--hpc-cap` added.
`index` and `search` must use the same `--hpc-cap`.

### 4. Calling

```bash
vespa call --reference ref.fa --bam tumor.ec.hp.bam --sfs tumor.dedup.sfs --normal-contigs-bam normal.ec.hp.bam \
  --bed-exclusion excluded.bed --clipped --no-clipped-fallback \
  --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --threads 16 > calls.unsorted.vcf
bcftools sort calls.unsorted.vcf -o calls.vcf
```
`--normal-contigs-bam` is the alignment of the corrected normal reads.
`--bed-exclusion` drops calls in the given regions (e.g., centromeres and telomeres).
`--min-cluster-weight` is the minimum number of tumor reads supporting a call: higher values increase precision at the cost of recall; we use 5 by default.

Calls that the normal reads show as germline are removed.
The other germline filters flag a call in the `FILTER` field instead of removing it: `GERMLINE`, `LOWPOWER` (too few normal reads to rule out a germline event), `POLYNORMAL` (the locus is length-polymorphic in the normal) and `REPEATDEPTH` (the normal depth is far above the median: a mismapping sink).
Somatic calls are those with `FILTER=PASS`.

## Reproducing the paper

The scripts that produce the callsets and the evaluations of the paper, with the download of their inputs, are in [VESPA-experiments](https://github.com/dcmonti/VESPA-experiments).
They use the release v1.0.1.

## License

MIT, see [LICENSE](LICENSE).
