#!/usr/bin/env python3
"""Deduplicate an SFS file and drop SFS longer than a cap.

SFS wire format (VESPA search): one line per SFS, tab-separated
    read_name  start  length  htag
where read_name is written only on the first SFS of a read and replaced by '*'
on the following ones. Duplicates are identified by (read, start, length).
The output keeps the same convention: the first SFS written for a read carries
its name, the others '*'. Header lines ('#...') are copied.

Usage: dedup_sfs.py in.sfs out.sfs [--max-length N]   (N = 0: no cap)
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sfs_in")
    ap.add_argument("sfs_out")
    ap.add_argument("--max-length", type=int, default=0, help="drop SFS longer than this (0 = no cap)")
    args = ap.parse_args()

    seen = set()
    written = set()
    read = ""
    n_in = n_dup = n_long = n_out = 0
    with open(args.sfs_in) as fin, open(args.sfs_out, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                fout.write(line + "\n")
                continue
            n_in += 1
            cols = line.split("\t")
            if cols[0] != "*":
                read = cols[0]
            key = (read, cols[1], cols[2])
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            if args.max_length and int(cols[2]) > args.max_length:
                n_long += 1
                continue
            if read in written:
                cols[0] = "*"
            else:
                cols[0] = read
                written.add(read)
            fout.write("\t".join(cols) + "\n")
            n_out += 1

    print(f"input SFS: {n_in}", file=sys.stderr)
    print(f"duplicates removed: {n_dup}", file=sys.stderr)
    print(f"longer than {args.max_length} bp removed: {n_long}", file=sys.stderr)
    print(f"output SFS: {n_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
