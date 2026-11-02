# ============================================================
# VESPA and hifiasm binaries (see tools.yaml and tools.mode in common.smk).
#
# Binaries are cached in tools.dir (default ~/.cache/vespa) under
# <tool>/<tag>/<mode>/, so they are fetched or built once and shared by runs.
#
# release: static x86-64 binaries attached to the GitHub release of the tag,
#          checked against the sha256 in tools.yaml.
# build:   the tag compiled from source with the system toolchain: git, make,
#          cmake >= 3.14, a C/C++ compiler with OpenMP, autotools, and the zlib,
#          bzip2 and xz development headers. VESPA fetches its other
#          dependencies (htslib, abPOA, ropebwt3, ...) at pinned commits.
# ============================================================


def _version_check_sh(binary, tool):
    """Shell snippet checking that `binary` reports the pinned tag. Returned through
    `params`, so braces are not re-formatted by Snakemake."""
    want = PINNED[tool]["tag"]
    strip = "sed 's/^VESPA //' | " if tool == "vespa" else ""
    on_mismatch = (
        'echo "WARNING: version mismatch allowed by tools.allow_version_mismatch" >&2'
        if TOOLS["allow_version_mismatch"] else "exit 1"
    )
    return (
        f'got=$({binary} --version 2>&1 | head -1 | {strip}awk "{{print \\$1}}")\n'
        f'if [ "$got" != "{want}" ]; then\n'
        f'    echo "{_EXE[tool]} reports version $got, tools.yaml pins {want}" >&2; {on_mismatch}\n'
        f'fi'
    )


def _release_url(tool):
    p = PINNED[tool]
    return f"{p['repo'].rstrip('/')}/releases/download/{p['tag']}/{p['asset']}"


if TOOLS_MODE == "release":

    for _tool, _bin in (("vespa", VESPA), ("hifiasm", HIFIASM)):

        rule:
            name:
                f"fetch_{_tool}"
            output:
                binary=_bin,
            params:
                url=_release_url(_tool),
                sha=PINNED[_tool].get("sha256") or "",
                check=_version_check_sh(_bin + ".tmp", _tool),
                tmp=_bin + ".tmp",
            localrule: True
            log:
                os.path.join(TOOLS_DIR, "logs", f"fetch_{_tool}.log"),
            shell:
                """
                mkdir -p $(dirname {output.binary})
                trap 'rm -f {params.tmp}' EXIT
                curl -fsSL --retry 3 -o {params.tmp} {params.url} 2> {log}
                if [ -n "{params.sha}" ]; then
                    echo "{params.sha}  {params.tmp}" | sha256sum -c - >> {log} 2>&1
                else
                    echo "WARNING: no sha256 in tools.yaml, download not verified" | tee -a {log} >&2
                fi
                chmod +x {params.tmp}
                {params.check}
                mv {params.tmp} {output.binary}
                """

elif TOOLS_MODE == "build":

    rule build_vespa:
        output:
            binary=VESPA,
        params:
            repo=PINNED["vespa"]["repo"],
            tag=PINNED["vespa"]["tag"],
            work=os.path.join(TOOLS_DIR, "vespa", PINNED["vespa"]["tag"], "build", "src"),
            check=_version_check_sh(VESPA, "vespa"),
        threads: THREADS["default"]
        log:
            os.path.join(TOOLS_DIR, "logs", "build_vespa.log"),
        shell:
            """
            rm -rf {params.work}
            git clone -q --depth 1 --branch {params.tag} {params.repo} {params.work} > {log} 2>&1
            cmake -S {params.work} -B {params.work}/build -DCMAKE_BUILD_TYPE=Release >> {log} 2>&1
            make -C {params.work}/build -j {threads} >> {log} 2>&1
            # VESPA writes its binary into the source directory
            cp {params.work}/vespa {output.binary}
            rm -rf {params.work}
            {params.check}
            """

    rule build_hifiasm:
        output:
            binary=HIFIASM,
        params:
            repo=PINNED["hifiasm"]["repo"],
            tag=PINNED["hifiasm"]["tag"],
            work=os.path.join(TOOLS_DIR, "hifiasm", PINNED["hifiasm"]["tag"], "build", "src"),
            check=_version_check_sh(HIFIASM, "hifiasm"),
        threads: THREADS["default"]
        log:
            os.path.join(TOOLS_DIR, "logs", "build_hifiasm.log"),
        shell:
            """
            rm -rf {params.work}
            git clone -q --depth 1 --branch {params.tag} {params.repo} {params.work} > {log} 2>&1
            make -C {params.work} -j {threads} >> {log} 2>&1
            cp {params.work}/hifiasm {output.binary}
            rm -rf {params.work}
            {params.check}
            """
