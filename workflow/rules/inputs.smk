# ============================================================
# Input normalization.
#
#   raw/{sample}.bam    raw reads aligned to any reference: needed only to split the
#                       error correction by chromosome (an aligned BAM input is used as is)
#   raw/{sample}.fq.gz  raw reads as FASTQ: those of `chromosomes` with split_by_chr,
#                       all of them otherwise. Input of whole-genome error correction
#                       and of the raw alignment used for haplotagging.
# ============================================================

LINK_SAMPLES = [s for s in SAMPLES if KIND[s] == "bam_aligned"]
ALIGN_SAMPLES = [s for s in SAMPLES if KIND[s] != "bam_aligned"]
BAM_SAMPLES = [s for s in SAMPLES if KIND[s].startswith("bam")]
READ_SAMPLES = [s for s in SAMPLES if not KIND[s].startswith("bam")]


def _one_of(samples):
    return "|".join(samples)


if SPLIT and LINK_SAMPLES:

    # Aligned BAM input: used as is, through a link, so its index can be created
    # in outdir when the input directory is read-only or has no index.
    rule link_raw_bam:
        input:
            bam=lambda w: INPUT[w.sample],
        output:
            bam=raw_bam("{sample}"),
            bai=raw_bam("{sample}") + ".bai",
        wildcard_constraints:
            sample=_one_of(LINK_SAMPLES),
        threads: THREADS["default"]
        log:
            log("link_raw_bam.{sample}"),
        conda:
            os.path.join(ENVS, "align.yml")
        shell:
            """
            ln -sf {input.bam} {output.bam}
            if [ -s {input.bam}.bai ]; then
                ln -sf {input.bam}.bai {output.bai}
            else
                samtools index -@ {threads} {output.bam} 2> {log}
            fi
            """


if SPLIT and ALIGN_SAMPLES:

    # Unaligned input: align the raw reads to the FULL reference (not the subset),
    # so that reads from alt contigs, decoys or other chromosomes land there and
    # are dropped by the per-chromosome split, instead of being forced onto
    # `chromosomes`.
    rule align_raw_for_split:
        input:
            reads=lambda w: INPUT[w.sample],
            ref=REF,
        output:
            bam=temp(raw_bam("{sample}")),
            bai=temp(raw_bam("{sample}") + ".bai"),
        wildcard_constraints:
            sample=_one_of(ALIGN_SAMPLES),
        params:
            # minimap2 reads FASTQ/FASTA (gzipped or not) directly; a uBAM is streamed
            reads=lambda w: "-" if KIND[w.sample] == "bam_unaligned" else INPUT[w.sample],
            ubam=lambda w: INPUT[w.sample] if KIND[w.sample] == "bam_unaligned" else "",
            tmp=lambda w: out("raw", f"{w.sample}.sort.tmp"),
        threads: THREADS["minimap2"]
        log:
            log("align_raw_for_split.{sample}"),
        benchmark:
            bench("align_raw_for_split.{sample}")
        conda:
            os.path.join(ENVS, "align.yml")
        shell:
            """
            if [ -n "{params.ubam}" ]; then samtools fastq -@ 4 {params.ubam} 2>/dev/null; fi \
              | minimap2 -ax map-hifi -t {threads} {input.ref} {params.reads} 2> {log} \
              | samtools sort -@ 8 -m 2G -T {params.tmp} -o {output.bam} - 2>> {log}
            samtools index -@ 8 {output.bam} 2>> {log}
            """


if not SPLIT and BAM_SAMPLES:

    rule bam_to_fastq:
        input:
            bam=lambda w: INPUT[w.sample],
        output:
            fq=temp(raw_reads("{sample}")),
        wildcard_constraints:
            sample=_one_of(BAM_SAMPLES),
        threads: THREADS["default"]
        log:
            log("bam_to_fastq.{sample}"),
        conda:
            os.path.join(ENVS, "align.yml")
        shell:
            """
            # -F 0x900: one record per read (no secondary / supplementary)
            samtools fastq -@ {threads} -F 0x900 {input.bam} 2> {log} | gzip -1 > {output.fq}
            """


def raw_fastq(s):
    """Raw reads given to hifiasm (whole genome) and to the raw alignment."""
    return INPUT[s] if not SPLIT and s in READ_SAMPLES else raw_reads(s)
