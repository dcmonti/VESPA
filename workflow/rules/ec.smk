# ============================================================
# Error correction with hifiasm.
#
#   tumor : --bin-only (only the corrected reads are needed)
#   normal: full run, because the r_utg graph provides the contigs of the index
#
# Thresholds: --het-cov for the normal; for the tumor either --het-cov or, with
# tumor_purity > 0, the negative-binomial model of the hifiasm fork (see its README).
#
# split_by_chr (default): reads are binned by the chromosome of their primary
# alignment and corrected one chromosome at a time. Below ~55x real depth
# whole-genome EC injects errors instead of removing them (chr20 at 27x: 1.97
# mismatches/kb vs 1.36 in the raw reads; 1.29 per chromosome). Reads that are
# unmapped or whose primary alignment is outside `chromosomes` are dropped.
#
# Outputs (both modes): ec/{sample}.ec.fa.gz and ec/normal.r_utg.fa.gz
# ============================================================

CHROM_RE = "|".join(re.escape(c) for c in CHROMS)


def _gfa_to_fasta(prefix):
    """awk program turning GFA segment lines into FASTA, names prefixed.
    Returned through `params`, so braces are NOT re-formatted by Snakemake."""
    return f"awk -v p='{prefix}' '/^S/{{print \">\" p $2; print $3}}'"


if SPLIT:

    rule ec_extract_chrom:
        input:
            bam=lambda w: raw_bam(w.sample),
            bai=lambda w: raw_bam(w.sample) + ".bai",
        output:
            fq=temp(out("ec", "chroms", "{sample}", "{chrom}.fq.gz")),
        wildcard_constraints:
            chrom=CHROM_RE,
        threads: 2
        log:
            log("ec_extract_chrom.{sample}.{chrom}"),
        conda:
            os.path.join(ENVS, "align.yml")
        shell:
            """
            # -F 0x900: each read once, on the chromosome of its primary alignment
            samtools view -u -F 0x900 {input.bam} {wildcards.chrom} 2> {log} \
              | samtools fastq - 2>> {log} | gzip -1 > {output.fq}
            """

    rule hifiasm_chrom:
        input:
            hifiasm_bin=HIFIASM_IN,
            fq=out("ec", "chroms", "{sample}", "{chrom}.fq.gz"),
        output:
            ec=temp(out("ec", "chroms", "{sample}", "{chrom}.ec.fa.gz")),
            contigs=temp(out("ec", "chroms", "{sample}", "{chrom}.r_utg.fa.gz")),
        wildcard_constraints:
            chrom=CHROM_RE,
        params:
            hifiasm=HIFIASM,
            flags=lambda w: ec_flags(w.sample),
            bin_only=lambda w: "--bin-only" if w.sample == "tumor" else "",
            extra=EC["extra"],
            work=lambda w: out("ec", "chroms", w.sample, f"{w.chrom}.work"),
            gfa2fa=lambda w: _gfa_to_fasta(f"{w.chrom}_"),
        threads: THREADS["hifiasm"]
        log:
            log("hifiasm.{sample}.{chrom}"),
        benchmark:
            bench("hifiasm.{sample}.{chrom}")
        shell:
            """
            rm -rf {params.work}; mkdir -p {params.work}
            if [ -z "$(gzip -dc {input.fq} | head -c 1)" ]; then
                # no read on this chromosome (e.g. chrY in a female sample)
                echo "no reads on {wildcards.chrom}" > {log}
                : | gzip > {output.ec}; : | gzip > {output.contigs}
            else
                {params.hifiasm} -t {threads} --write-ec {params.flags} {params.extra} \
                    -o {params.work}/ec {params.bin_only} {input.fq} > {log} 2>&1
                [ -s {params.work}/ec.ec.fa ] || {{ echo "hifiasm wrote no ec.fa" >> {log}; exit 1; }}
                gzip -1 -c {params.work}/ec.ec.fa > {output.ec}
                # hifiasm reuses utg000001l in every run: prefix contigs with the chromosome
                gfa={params.work}/ec.bp.r_utg.gfa
                if [ -s "$gfa" ]; then {params.gfa2fa} "$gfa" | gzip -1 > {output.contigs}
                else : | gzip > {output.contigs}; fi
            fi
            rm -rf {params.work}
            """

    rule ec_gather:
        input:
            ec=lambda w: expand(out("ec", "chroms", w.sample, "{chrom}.ec.fa.gz"), chrom=CHROMS),
        output:
            ec=out("ec", "{sample}.ec.fa.gz"),
        wildcard_constraints:
            sample="tumor|normal",
        shell:
            # concatenated gzip members are a valid gzip stream (read by minimap2 and vespa)
            "cat {input.ec} > {output.ec}"

    rule raw_gather:
        input:
            fq=lambda w: expand(out("ec", "chroms", w.sample, "{chrom}.fq.gz"), chrom=CHROMS),
        output:
            fq=temp(raw_reads("{sample}")),
        wildcard_constraints:
            sample="tumor|normal",
        shell:
            # raw reads of `chromosomes`, for the raw alignment used by haplotagging
            "cat {input.fq} > {output.fq}"

    rule contigs_gather:
        input:
            contigs=expand(out("ec", "chroms", "normal", "{chrom}.r_utg.fa.gz"), chrom=CHROMS),
        output:
            contigs=NORMAL_CONTIGS,
        shell:
            "cat {input.contigs} > {output.contigs}"

else:

    rule hifiasm_wgs:
        input:
            hifiasm_bin=HIFIASM_IN,
            reads=lambda w: raw_fastq(w.sample),
        output:
            ec=out("ec", "{sample}.ec.fa.gz"),
            contigs=out("ec", "{sample}.r_utg.fa.gz"),
        wildcard_constraints:
            sample="tumor|normal",
        params:
            hifiasm=HIFIASM,
            flags=lambda w: ec_flags(w.sample),
            bin_only=lambda w: "--bin-only" if w.sample == "tumor" else "",
            extra=EC["extra"],
            work=lambda w: out("ec", f"{w.sample}.work"),
            gfa2fa=_gfa_to_fasta(""),
        threads: THREADS["hifiasm"]
        log:
            log("hifiasm.{sample}"),
        benchmark:
            bench("hifiasm.{sample}")
        shell:
            """
            rm -rf {params.work}; mkdir -p {params.work}
            {params.hifiasm} -t {threads} --write-ec {params.flags} {params.extra} \
                -o {params.work}/ec {params.bin_only} {input.reads} > {log} 2>&1
            [ -s {params.work}/ec.ec.fa ] || {{ echo "hifiasm wrote no ec.fa" >> {log}; exit 1; }}
            gzip -1 -c {params.work}/ec.ec.fa > {output.ec}
            gfa={params.work}/ec.bp.r_utg.gfa
            if [ -s "$gfa" ]; then {params.gfa2fa} "$gfa" | gzip -1 > {output.contigs}
            else : | gzip > {output.contigs}; fi
            rm -rf {params.work}
            """
