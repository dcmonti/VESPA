# ============================================================
# FMD index, built in two phases so that it absorbs the germline sequence of the
# normal:
#   1. base  = reference (all its contigs) + normal r_utg contigs
#   2. search the normal EC reads against the base index: the reads that still
#      produce SFS carry sequence the contigs miss ("missing" reads)
#   3. full  = base + missing reads, appended with ropebwt3's -i. A BWT is
#      canonical for its sequence multiset, so appending is equivalent to a
#      rebuild, without re-indexing the reference.
# -d (fermi-delta format) is mandatory: `vespa search` only reads FMD.
# With index.hpc every step homopolymer-compresses its input and vespa writes a
# `.hpc` sidecar next to the index.
# ============================================================


def _index_checks(fmd, logf):
    """Shell snippet failing when ropebwt3 did not complete. (Its "inserted N symbols"
    line only appears when batches are merged, so it cannot be required.)"""
    hpc = f'[ -f {fmd}.hpc ] || {{ echo "missing {fmd}.hpc sidecar" >&2; exit 1; }}' if IDX["hpc"] else ":"
    return (
        f'grep -qa "All done!" {logf} || {{ echo "index build did not complete, see {logf}" >&2; exit 1; }}\n'
        f'[ -s {fmd} ] || {{ echo "{fmd} is empty, see {logf}" >&2; exit 1; }}\n'
        + hpc
    )


def _search_checks(sfs, logf):
    return (
        f'grep -qa "All done!" {logf} || {{ echo "search did not complete, see {logf}" >&2; exit 1; }}\n'
        f'[ "$(grep -vc "^#" {sfs} || true)" -gt 0 ] || {{ echo "{sfs} has no SFS, see {logf}" >&2; exit 1; }}'
    )


rule index_base:
    input:
        vespa_bin=VESPA_IN,
        ref=REF,
        contigs=NORMAL_CONTIGS,
    output:
        fmd_outputs(BASE_FMD),
    params:
        vespa=VESPA,
        batch=IDX["batch_size"],
        sais=IDX["sais_threads"],
        hpc=HPC_FLAG,
        checks=_index_checks(BASE_FMD, log("index_base")),
    threads: THREADS["vespa"]
    log:
        log("index_base"),
    benchmark:
        bench("index_base")
    shell:
        """
        {params.vespa} index -m {params.batch} -t {threads} -p {params.sais} -d {params.hpc} \
            {input.ref} {input.contigs} -o {output[0]} > {log} 2>&1
        {params.checks}
        """


rule search_normal:
    input:
        vespa_bin=VESPA_IN,
        fmd=fmd_outputs(BASE_FMD),
        reads=ec_reads("normal"),
    output:
        sfs=NORMAL_SFS,
    params:
        vespa=VESPA,
        bsize=SRCH["bsize"],
        hpc=HPC_FLAG,
        checks=_search_checks(NORMAL_SFS, log("search_normal")),
    threads: THREADS["vespa"]
    log:
        log("search_normal"),
    benchmark:
        bench("search_normal")
    shell:
        """
        {params.vespa} search --index {input.fmd[0]} --fastx {input.reads} \
            --bsize {params.bsize} --threads {threads} {params.hpc} > {output.sfs} 2> {log}
        {params.checks}
        """


rule extract_missing_reads:
    input:
        sfs=NORMAL_SFS,
        reads=ec_reads("normal"),
    output:
        fa=MISSING_READS,
    log:
        log("extract_missing_reads"),
    conda:
        os.path.join(ENVS, "python.yml")
    shell:
        "python3 {SCRIPTS}/extract_missing_reads.py {input.sfs} {input.reads} {output.fa} 2> {log}"


rule index_full:
    input:
        vespa_bin=VESPA_IN,
        base=fmd_outputs(BASE_FMD),
        missing=MISSING_READS,
    output:
        fmd_outputs(FULL_FMD),
    params:
        vespa=VESPA,
        batch=IDX["batch_size"],
        sais=IDX["sais_threads"],
        hpc=HPC_FLAG,
        checks=_index_checks(FULL_FMD, log("index_full")),
    threads: THREADS["vespa"]
    log:
        log("index_full"),
    benchmark:
        bench("index_full")
    shell:
        """
        {params.vespa} index -m {params.batch} -t {threads} -p {params.sais} -d \
            -i {input.base[0]} {params.hpc} {input.missing} -o {output[0]} > {log} 2>&1
        {params.checks}
        """
