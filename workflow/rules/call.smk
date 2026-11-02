# ============================================================
# vespa call.
#
# The normal EC alignments feed the in-call germline filter and the post-call
# gates (FILTER GERMLINE, LOWPOWER, POLYNORMAL, REPEATDEPTH). HP tags, when
# present, are read from the tumor BAM for the per-haplotype split before POA,
# and from both BAMs by the gates.
# ============================================================


rule call:
    input:
        vespa_bin=VESPA_IN,
        ref=CALL_REF,
        tumor=call_bam("tumor"),
        normal=call_bam("normal"),
        sfs=TUMOR_SFS_FINAL,
    output:
        vcf=temp(out("call", "calls.unsorted.vcf")),
        **({"clusters": out("call", "calls.clusters")} if CALL["save_clusters"] else {}),
        **(
            {"clip_clusters": out("call", "calls.clipped_clusters")}
            if CALL["save_clusters"] and CALL["clipped"] else {}
        ),
    params:
        vespa=VESPA,
        weight=CALL["min_cluster_weight"],
        min_len=CALL["min_sv_length"],
        mapq=CALL["min_mapq"],
        clipped="--clipped --no-clipped-fallback" if CALL["clipped"] else "",
        noht="--noht" if CALL["noht"] else "",
        bed=f"--bed-exclusion {EXCL_BED}" if EXCL_BED else "",
        dumps=(
            (f"--clusters {out('call', 'calls.clusters')}" if CALL["save_clusters"] else "")
            + (f" --clip-clusters {out('call', 'calls.clipped_clusters')}"
               if CALL["save_clusters"] and CALL["clipped"] else "")
        ),
        extra=CALL["extra"],
    threads: THREADS["vespa"]
    log:
        log("call"),
    benchmark:
        bench("call")
    shell:
        """
        {params.vespa} call --reference {input.ref} --bam {input.tumor} --sfs {input.sfs} \
            --normal-contigs-bam {input.normal} {params.bed} \
            --min-cluster-weight {params.weight} --min-sv-length {params.min_len} \
            --min-mapq {params.mapq} {params.clipped} {params.noht} {params.dumps} \
            {params.extra} --threads {threads} > {output.vcf} 2> {log}
        """


rule sort_vcf:
    input:
        vcf=out("call", "calls.unsorted.vcf"),
    output:
        vcf=CALLS_VCF,
        pass_vcf=CALLS_PASS,
    log:
        log("sort_vcf"),
    conda:
        os.path.join(ENVS, "whatshap.yml")
    shell:
        """
        bcftools sort {input.vcf} -o {output.vcf} 2> {log}
        bcftools view -f PASS {output.vcf} -o {output.pass_vcf} 2>> {log}
        """
