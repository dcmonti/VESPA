# ============================================================
# Alignments to the call reference, minimap2 HiFi preset at its defaults:
#   align/{sample}.ec.bam   error-corrected reads, given to vespa call
#   align/{sample}.raw.bam  raw reads, used only for haplotagging (DeepVariant needs
#                           base qualities, which the error-corrected reads lack)
# ============================================================


def _align_shell(reads):
    return (
        f"minimap2 -ax map-hifi -t {{threads}} {{input.ref}} {reads} 2> {{log}} \\\n"
        "  | samtools sort -@ 8 -m 4G -T {params.tmp} -o {output.bam} - 2>> {log}\n"
        "samtools index -@ 8 {output.bam} 2>> {log}\n"
    )


rule align_ec:
    input:
        reads=lambda w: ec_reads(w.sample),
        ref=CALL_REF,
    output:
        # with haplotagging only the tagged copy (align/{sample}.ec.hp.bam) is kept
        bam=temp(ec_bam("{sample}")) if HT_MODE != "none" else ec_bam("{sample}"),
        bai=temp(ec_bam("{sample}") + ".bai") if HT_MODE != "none" else ec_bam("{sample}") + ".bai",
    wildcard_constraints:
        sample="tumor|normal",
    params:
        tmp=lambda w: out("align", f"{w.sample}.ec.sort.tmp"),
    threads: THREADS["minimap2"]
    log:
        log("align_ec.{sample}"),
    benchmark:
        bench("align_ec.{sample}")
    conda:
        os.path.join(ENVS, "align.yml")
    shell:
        _align_shell("{input.reads}")


if HT_MODE != "none":

    rule align_raw:
        input:
            reads=lambda w: raw_fastq(w.sample),
            ref=CALL_REF,
        output:
            bam=temp(raw_call_bam("{sample}")),
            bai=temp(raw_call_bam("{sample}") + ".bai"),
        wildcard_constraints:
            sample="tumor|normal",
        params:
            tmp=lambda w: out("align", f"{w.sample}.raw.sort.tmp"),
        threads: THREADS["minimap2"]
        log:
            log("align_raw.{sample}"),
        benchmark:
            bench("align_raw.{sample}")
        conda:
            os.path.join(ENVS, "align.yml")
        shell:
            _align_shell("{input.reads}")
