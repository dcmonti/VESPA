# ============================================================
# Configuration: defaults, validation, input detection, paths.
# Everything here runs at DAG-construction time.
# ============================================================

import copy
import os
import re
import shutil
import sys
import subprocess

from snakemake.exceptions import WorkflowError

# ---- Defaults: every key the workflow reads is listed here ----
DEFAULTS = {
    "tools": {
        "mode": "release",
        "dir": "",
        "vespa_bin": "",
        "hifiasm_bin": "",
        "allow_version_mismatch": False,
    },
    "tumor": "",
    "normal": "",
    "reference": "",
    "chromosomes": [],
    "exclusion_bed": "",
    "outdir": "results",
    "error_correction": {
        "split_by_chr": True,
        "normal_het_cov": 4,
        "tumor_het_cov": 6,
        "tumor_purity": 0,
        "tumor_clonality": 1.0,
        "tumor_dispersion": 0.02,
        "tumor_somatic_fnr": 0.10,
        "extra": "",
    },
    "haplotag": {
        "mode": "deepvariant",
        "phased_vcf": "",
        "container_runtime": "singularity",
        "deepvariant_image": "docker://google/deepvariant:1.9.0",
        "min_tagged_fraction": 0.2,
    },
    "index": {
        "hpc": True,
        "hpc_cap": 2,
        "batch_size": "4G",
        "sais_threads": 0,
    },
    "search": {
        "bsize": 10000,
        "max_sfs_length": 500,
    },
    "call": {
        "min_cluster_weight": 5,
        "min_sv_length": 50,
        "min_mapq": 60,
        "clipped": True,
        "noht": False,
        "save_clusters": False,
        "extra": "",
    },
    "threads": {
        "default": 8,
        "hifiasm": 32,
        "minimap2": 32,
        "deepvariant": 16,
        "whatshap": 4,
        "vespa": 32,
    },
}


def _merge(defaults, user, path=""):
    """Recursively overlay `user` on `defaults`, rejecting unknown keys."""
    out = copy.deepcopy(defaults)
    for key, value in user.items():
        where = f"{path}{key}"
        if key not in defaults:
            raise WorkflowError(f"Unknown config key '{where}'.")
        if isinstance(defaults[key], dict):
            if not isinstance(value, dict):
                raise WorkflowError(f"Config key '{where}' must be a mapping.")
            out[key] = _merge(defaults[key], value, where + ".")
        else:
            out[key] = value
    return out


CFG = _merge(DEFAULTS, dict(config))
EC = CFG["error_correction"]
HT = CFG["haplotag"]
IDX = CFG["index"]
SRCH = CFG["search"]
CALL = CFG["call"]
THREADS = CFG["threads"]

SAMPLES = ["tumor", "normal"]
OUT = CFG["outdir"]


def out(*parts):
    return os.path.join(OUT, *parts)


def _die(msg):
    raise WorkflowError(msg)


def _run(cmd):
    """Run a command at DAG time and return stdout+stderr (never raises)."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return p.stdout + p.stderr
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"__failed__: {e}"


# ---- Binaries ----
# The versions the workflow runs are pinned in tools.yaml. tools.mode:
#   release  download the static binaries of the pinned release (default)
#   build    compile the pinned tags from source
#   path     use tools.vespa_bin / tools.hifiasm_bin as given (development)
# In every mode `--version` must report the pinned tag.
import yaml

with open(os.path.join(workflow.basedir, "tools.yaml")) as _f:
    PINNED = yaml.safe_load(_f)

TOOLS = CFG["tools"]
TOOLS_MODE = TOOLS["mode"]
if TOOLS_MODE not in ("release", "build", "path"):
    _die("tools.mode must be release, build or path.")
TOOLS_DIR = os.path.abspath(os.path.expanduser(TOOLS["dir"] or "~/.cache/vespa"))
_EXE = {"vespa": "vespa", "hifiasm": "hifiasm"}


def version_token(binary, tool):
    """Version in `--version` ('VESPA v1.0.1 (...)' -> 'v1.0.1', 'v1.0.0-somatic (...)' -> 'v1.0.0-somatic')."""
    text = _run([binary, "--version"]).strip().splitlines()
    line = text[0] if text else ""
    if tool == "vespa" and line.startswith("VESPA "):
        line = line[len("VESPA "):]
    parts = line.split()
    return parts[0] if parts else ""


def check_version(binary, tool):
    """Stop (or warn) when `binary` does not report the pinned tag."""
    want = PINNED[tool]["tag"]
    got = version_token(binary, tool)
    if got == want:
        return
    msg = (f"{_EXE[tool]} '{binary}' reports version '{got}', the workflow pins '{want}' "
           f"(tools.yaml).")
    if TOOLS["allow_version_mismatch"]:
        print(f"WARNING: {msg} Continuing because tools.allow_version_mismatch is set.",
              file=sys.stderr)
    else:
        _die(msg + " Use tools.mode: release/build, or set tools.allow_version_mismatch "
                   "for a development binary.")


def tool_path(tool):
    if TOOLS_MODE == "path":
        b = TOOLS[f"{tool}_bin"]
        if not b:
            _die(f"tools.mode is 'path' but tools.{tool}_bin is empty.")
        b = shutil.which(b) or b
        if not os.access(b, os.X_OK):
            _die(f"tools.{tool}_bin '{b}' not found or not executable.")
        return os.path.abspath(b)
    return os.path.join(TOOLS_DIR, tool, PINNED[tool]["tag"], TOOLS_MODE, _EXE[tool])


VESPA = tool_path("vespa")
HIFIASM = tool_path("hifiasm")
for _tool, _bin in (("vespa", VESPA), ("hifiasm", HIFIASM)):
    # release/build binaries are checked by the rule that produces them; a cached
    # one is re-checked here, so a stale cache cannot go unnoticed.
    if TOOLS_MODE == "path" or os.path.exists(_bin):
        check_version(_bin, _tool)

# As rule inputs: a user binary is `ancient`, so rebuilding it does not rerun the
# whole workflow.
VESPA_IN = ancient(VESPA) if TOOLS_MODE == "path" else VESPA
HIFIASM_IN = ancient(HIFIASM) if TOOLS_MODE == "path" else HIFIASM

# ---- Inputs ----
FASTQ_EXT = (".fq", ".fastq", ".fq.gz", ".fastq.gz")
FASTA_EXT = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def detect_input(path):
    """Return 'bam_aligned', 'bam_unaligned', 'fastq' or 'fasta'."""
    if not path:
        _die("Both 'tumor' and 'normal' must be set in the config.")
    if not isinstance(path, str):
        _die(f"Input '{path}': give a single file per sample.")
    if not os.path.exists(path):
        _die(f"Input file '{path}' not found.")
    low = path.lower()
    if low.endswith(".cram"):
        _die(f"Input '{path}': CRAM is not supported, convert it to BAM first.")
    if low.endswith(FASTQ_EXT):
        return "fastq"
    if low.endswith(FASTA_EXT):
        return "fasta"
    if low.endswith(".bam"):
        header = _run(["samtools", "view", "-H", path])
        if "@SQ" not in header:
            return "bam_unaligned"
        # One mapped record is enough; `-c` would scan the whole file.
        p = subprocess.Popen(["samtools", "view", "-F", "4", path],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        first = p.stdout.readline()
        p.kill()
        return "bam_aligned" if first else "bam_unaligned"
    _die(f"Input '{path}': unrecognised extension (expected .bam, FASTQ or FASTA).")


INPUT = {s: os.path.abspath(CFG[s]) if CFG[s] else "" for s in SAMPLES}
KIND = {s: detect_input(CFG[s]) for s in SAMPLES}


def bam_contigs(path):
    return [l.split("\t")[1][3:] for l in _run(["samtools", "view", "-H", path]).splitlines()
            if l.startswith("@SQ")]


# ---- Reference and chromosomes ----
REF = CFG["reference"]
if not REF or not os.path.exists(REF):
    _die(f"reference '{REF}' not found.")
if not os.path.exists(REF + ".fai"):
    if "__failed__" in _run(["samtools", "faidx", REF]) or not os.path.exists(REF + ".fai"):
        _die(f"Cannot index '{REF}': use an uncompressed or bgzip-compressed FASTA.")
REF_CONTIGS = [l.split("\t")[0] for l in open(REF + ".fai")]

PRIMARY_RE = re.compile(r"^(chr)?([0-9]+|X|Y)$")
if CFG["chromosomes"]:
    CHROMS = [str(c) for c in CFG["chromosomes"]]
    _missing = [c for c in CHROMS if c not in REF_CONTIGS]
    if _missing:
        _die(f"chromosomes not in {REF}.fai: {', '.join(_missing)}")
else:
    CHROMS = [c for c in REF_CONTIGS if PRIMARY_RE.match(c)]
    if not CHROMS:
        _die(f"No primary chromosome found in {REF}.fai: set 'chromosomes' explicitly.")

SPLIT = bool(EC["split_by_chr"])
for s in SAMPLES:
    if SPLIT and KIND[s] == "bam_aligned":
        _missing = [c for c in CHROMS if c not in set(bam_contigs(INPUT[s]))]
        if _missing:
            _more = " ..." if len(_missing) > 5 else ""
            _die(
                f"{s} BAM '{INPUT[s]}' has no {', '.join(_missing[:5])}"
                f"{_more}: chromosome names must match "
                "'chromosomes' (e.g. chr1 vs 1). The BAM does not need to be aligned to "
                "'reference', it is only used to bin reads by chromosome."
            )

# Reads are aligned, haplotagged and called on `reference` restricted to `chromosomes`, so
# that they cannot land on alt contigs or decoys and lose mapping quality. The index is
# built from the whole `reference`.
CALL_REF = out("ref", "call.fa")

EXCL_BED = CFG["exclusion_bed"]
if EXCL_BED and not os.path.exists(EXCL_BED):
    _die(f"exclusion_bed '{EXCL_BED}' not found.")

# ---- Haplotagging ----
HT_MODE = HT["mode"]
if HT_MODE not in ("none", "vcf", "deepvariant"):
    _die("haplotag.mode must be none, vcf or deepvariant.")
if HT_MODE == "vcf":
    if not HT["phased_vcf"] or not os.path.exists(HT["phased_vcf"]):
        _die("haplotag.mode is 'vcf' but haplotag.phased_vcf is missing.")
if HT_MODE == "deepvariant" and HT["container_runtime"] not in ("singularity", "apptainer", "docker"):
    _die("haplotag.container_runtime must be singularity, apptainer or docker.")

# ---- HPC flags: one definition shared by index and search ----
HPC_FLAG = f"--hpc --hpc-cap {int(IDX['hpc_cap'])}" if IDX["hpc"] else ""


# ---- EC threshold flags ----
def tumor_threshold_flags():
    if EC["tumor_purity"] and float(EC["tumor_purity"]) > 0:
        return (
            f"--tumor-purity {EC['tumor_purity']} --clonality {EC['tumor_clonality']} "
            f"--dispersion {EC['tumor_dispersion']} --somatic-fnr {EC['tumor_somatic_fnr']}"
        )
    return f"--het-cov {EC['tumor_het_cov']}"


def ec_flags(sample):
    return f"--het-cov {EC['normal_het_cov']}" if sample == "normal" else tumor_threshold_flags()


# ---- Paths ----
def raw_bam(s):
    """Raw reads as given, aligned to any reference: only to bin reads by chromosome."""
    return out("raw", f"{s}.bam")


def raw_reads(s):
    """Raw reads of `chromosomes` (split_by_chr) or all raw reads, as FASTQ."""
    return out("raw", f"{s}.fq.gz")


def raw_call_bam(s):
    """Raw reads aligned to CALL_REF: input of DeepVariant and whatshap."""
    return out("align", f"{s}.raw.bam")


def ec_reads(s):
    return out("ec", f"{s}.ec.fa.gz")


NORMAL_CONTIGS = out("ec", "normal.r_utg.fa.gz")


def ec_bam(s):
    return out("align", f"{s}.ec.bam")


def call_bam(s):
    """BAM given to vespa call: haplotagged when haplotagging is on."""
    return out("align", f"{s}.ec.hp.bam") if HT_MODE != "none" else ec_bam(s)


BASE_FMD = out("index", "base.fmd")
FULL_FMD = out("index", "full.fmd")
NORMAL_SFS = out("search", "normal.sfs")
MISSING_READS = out("search", "normal.missing.fa.gz")
TUMOR_SFS = out("search", "tumor.sfs")
TUMOR_SFS_FINAL = out("search", "tumor.final.sfs")
CALLS_VCF = out("call", "calls.vcf")
CALLS_PASS = out("call", "calls.PASS.vcf")


def fmd_outputs(fmd):
    return [fmd, fmd + ".hpc"] if IDX["hpc"] else [fmd]


def log(name):
    return out("logs", f"{name}.log")


def bench(name):
    return out("benchmarks", f"{name}.tsv")


SCRIPTS = os.path.join(workflow.basedir, "..", "scripts")
ENVS = os.path.join(workflow.basedir, "envs")

print(
    f"[vespa] tumor={KIND['tumor']} normal={KIND['normal']} "
    f"split_by_chr={SPLIT} haplotag={HT_MODE} "
    f"chromosomes={len(CHROMS)} ({CHROMS[0]}..{CHROMS[-1]})",
    file=sys.stderr,
)
