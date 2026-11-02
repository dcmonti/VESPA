# ============================================================
# Tumor search: sample-specific strings (SFS) of the tumor EC reads against the
# full index, then deduplication and length cap.
# ============================================================


rule search_tumor:
    input:
        vespa_bin=VESPA_IN,
        fmd=fmd_outputs(FULL_FMD),
        reads=ec_reads("tumor"),
    output:
        sfs=TUMOR_SFS,
    params:
        vespa=VESPA,
        bsize=SRCH["bsize"],
        hpc=HPC_FLAG,
        checks=_search_checks(TUMOR_SFS, log("search_tumor")),
    threads: THREADS["vespa"]
    log:
        log("search_tumor"),
    benchmark:
        bench("search_tumor")
    shell:
        """
        {params.vespa} search --index {input.fmd[0]} --fastx {input.reads} \
            --bsize {params.bsize} --threads {threads} {params.hpc} > {output.sfs} 2> {log}
        {params.checks}
        """


rule dedup_sfs:
    input:
        sfs=TUMOR_SFS,
    output:
        sfs=TUMOR_SFS_FINAL,
    params:
        max_length=SRCH["max_sfs_length"],
    log:
        log("dedup_sfs"),
    conda:
        os.path.join(ENVS, "python.yml")
    shell:
        "python3 {SCRIPTS}/dedup_sfs.py {input.sfs} {output.sfs} --max-length {params.max_length} 2> {log}"
