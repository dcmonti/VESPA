#!/usr/bin/env python3
"""Extract the normal EC reads that still produce SFS against the base index.

These reads carry sequence the r_utg contigs do not represent; appending them to
the index lets it absorb that germline sequence.

Usage: extract_missing_reads.py normal.sfs normal.ec.fa[.gz] missing.fa.gz
"""
import gzip
import sys


def main():
    sfs, reads, out = sys.argv[1:4]
    wanted = set()
    with open(sfs) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            name = line.split("\t", 1)[0]
            if name != "*":  # '*' repeats the previous read name
                wanted.add(name)

    opener = gzip.open if reads.endswith(".gz") else open
    n = 0
    keep = False
    with opener(reads, "rt") as fin, gzip.open(out, "wt", compresslevel=1) as fout:
        for line in fin:
            if line.startswith(">"):
                keep = line[1:].split(None, 1)[0] in wanted
                n += keep
            if keep:
                fout.write(line)

    print(f"read names with SFS: {len(wanted)}", file=sys.stderr)
    print(f"reads extracted: {n}", file=sys.stderr)
    if n == 0:
        sys.exit("no missing read extracted: check that the SFS and the reads come from the same run")


if __name__ == "__main__":
    main()
