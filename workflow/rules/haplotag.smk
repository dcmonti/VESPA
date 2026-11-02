# ============================================================
# Haplotagging (default: mode deepvariant).
#
#   mode deepvariant: DeepVariant on the raw normal alignment -> whatshap phase
#   mode vcf        : a phased VCF supplied by the user
# then, for both samples: whatshap haplotag on the raw alignment -> read list ->
# HP and PS tags copied onto the error-corrected alignment by read name (an EC read
# keeps the name of its raw read): align/{sample}.ec.hp.bam.
# ============================================================

if HT_MODE != "none":

    def ht_bam(s):
        return raw_call_bam(s)

    HT_REF = CALL_REF
    PHASED_VCF = out("haplotag", "normal.phased.vcf.gz")

    if HT_MODE == "deepvariant":

        def _container_cmd():
            binds = sorted({os.path.dirname(os.path.abspath(p)) for p in
                            [HT_REF, ht_bam("normal"), OUT + "/haplotag/x"]} | {os.path.abspath(OUT)})
            image = HT["deepvariant_image"]
            if HT["container_runtime"] == "docker":
                vols = " ".join(f"-v {b}:{b}" for b in binds)
                return f"docker run --rm {vols} -e TMPDIR={os.path.abspath(OUT)}/haplotag/tmp " \
                       f"{image.replace('docker://', '')}"
            vols = " ".join(f"-B {b}" for b in binds)
            return f"{HT['container_runtime']} exec --cleanenv {vols} " \
                   f"--env TMPDIR={os.path.abspath(OUT)}/haplotag/tmp {image}"

        rule deepvariant:
            input:
                bam=ht_bam("normal"),
                bai=ht_bam("normal") + ".bai",
                ref=HT_REF,
            output:
                vcf=out("haplotag", "normal.dv.vcf.gz"),
            params:
                container=_container_cmd(),
                inter=lambda w: os.path.abspath(out("haplotag", "dv_intermediate")),
                tmp=lambda w: os.path.abspath(out("haplotag", "tmp")),
                ref=os.path.abspath(HT_REF),
                bam=os.path.abspath(ht_bam("normal")),
                vcf=os.path.abspath(out("haplotag", "normal.dv.vcf.gz")),
            threads: THREADS["deepvariant"]
            log:
                log("deepvariant"),
            benchmark:
                bench("deepvariant")
            shell:
                """
                mkdir -p {params.inter} {params.tmp}
                {params.container} run_deepvariant --model_type=PACBIO \
                    --ref={params.ref} --reads={params.bam} \
                    --output_vcf={params.vcf} --num_shards={threads} \
                    --intermediate_results_dir={params.inter} > {log} 2>&1
                rm -rf {params.inter} {params.tmp}
                """

        rule whatshap_phase:
            input:
                vcf=out("haplotag", "normal.dv.vcf.gz"),
                bam=ht_bam("normal"),
                ref=HT_REF,
            output:
                vcf=PHASED_VCF,
                tbi=PHASED_VCF + ".tbi",
            log:
                log("whatshap_phase"),
            benchmark:
                bench("whatshap_phase")
            conda:
                os.path.join(ENVS, "whatshap.yml")
            shell:
                """
                whatshap phase -o {output.vcf} --reference {input.ref} --ignore-read-groups \
                    {input.vcf} {input.bam} > {log} 2>&1
                tabix -f -p vcf {output.vcf}
                """

    else:

        rule link_phased_vcf:
            input:
                vcf=HT["phased_vcf"],
            output:
                vcf=PHASED_VCF,
                tbi=PHASED_VCF + ".tbi",
            params:
                vcf=os.path.abspath(HT["phased_vcf"]),
            conda:
                os.path.join(ENVS, "whatshap.yml")
            shell:
                """
                ln -sf {params.vcf} {output.vcf}
                if [ -s {params.vcf}.tbi ]; then ln -sf {params.vcf}.tbi {output.tbi}
                else tabix -f -p vcf {output.vcf}; fi
                """

    rule whatshap_haplotag:
        input:
            vcf=PHASED_VCF,
            tbi=PHASED_VCF + ".tbi",
            bam=lambda w: ht_bam(w.sample),
            bai=lambda w: ht_bam(w.sample) + ".bai",
            ref=HT_REF,
        output:
            tsv=out("haplotag", "{sample}.haplotags.tsv"),
        wildcard_constraints:
            sample="tumor|normal",
        params:
            min_frac=HT["min_tagged_fraction"],
        threads: THREADS["whatshap"]
        log:
            log("whatshap_haplotag.{sample}"),
        benchmark:
            bench("whatshap_haplotag.{sample}")
        conda:
            os.path.join(ENVS, "whatshap.yml")
        shell:
            """
            # --tag-supplementary: supplementary alignments carry the SV evidence
            # --skip-missing-contigs: contigs absent from the VCF are expected
            whatshap haplotag --reference {input.ref} --output /dev/null \
                --output-haplotag-list {output.tsv} --ignore-read-groups \
                --tag-supplementary --skip-missing-contigs --output-threads={threads} \
                {input.vcf} {input.bam} > {log} 2>&1
            # whatshap tags almost nothing, silently, when contig names do not match
            awk -F'\\t' -v m={params.min_frac} '!/^#/ {{n++; if ($2=="H1"||$2=="H2") t++}}
                END {{ f = n ? t/n : 0; printf "tagged %d of %d reads (%.3f)\\n", t, n, f > "/dev/stderr";
                       if (f < m) {{ print "tagged fraction below haplotag.min_tagged_fraction: check that the VCF and the BAM use the same contigs" > "/dev/stderr"; exit 1 }} }}' \
                {output.tsv} 2>> {log}
            """

    rule tag_bam:
        input:
            bam=lambda w: ec_bam(w.sample),
            tsv=out("haplotag", "{sample}.haplotags.tsv"),
        output:
            bam=out("align", "{sample}.ec.hp.bam"),
            bai=out("align", "{sample}.ec.hp.bam.bai"),
        wildcard_constraints:
            sample="tumor|normal",
        threads: THREADS["default"]
        log:
            log("tag_bam.{sample}"),
        conda:
            os.path.join(ENVS, "python.yml")
        shell:
            """
            python3 {SCRIPTS}/tag_bam_from_list.py {input.bam} {input.tsv} {output.bam} \
                --threads {threads} --with-ps 2> {log}
            samtools index -@ {threads} {output.bam} 2>> {log}
            """
