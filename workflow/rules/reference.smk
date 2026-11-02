# ============================================================
# Reference for alignment, haplotagging and calling: `reference` restricted to
# `chromosomes`, so that reads corrected from those chromosomes cannot land on alt
# contigs, decoys or chrY and lose MAPQ (call runs with --min-mapq 60 by default).
# ============================================================


rule reference_subset:
    input:
        fa=REF,
    output:
        fa=CALL_REF,
        fai=CALL_REF + ".fai",
    params:
        chroms=" ".join(CHROMS),
    log:
        log("reference_subset"),
    conda:
        os.path.join(ENVS, "align.yml")
    shell:
        """
        samtools faidx {input.fa} {params.chroms} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """
