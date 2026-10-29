#!/usr/bin/env python3
"""Apply HP tags to a BAM from a whatshap haplotag read list.

The list is the one written by `whatshap haplotag --output-haplotag-list`
(header '#readname', columns: readname, haplotype H1/H2/none, phaseset,
chromosome). Matching is by read name, so the BAM to tag may hold different
alignments of the same reads (e.g. error-corrected reads).

Existing HP/PS tags are removed first: a read absent from the list, or listed
as 'none', ends up untagged instead of keeping a stale tag.

Only HP is needed downstream (VESPA reads HP); PS is optional (--with-ps).

Usage: tag_bam_from_list.py in.bam haplotags.tsv out.bam [--threads N] [--with-ps]
"""
import argparse
import sys

import pysam


def load_haplotags(tsv, with_ps):
    """read name -> HP (1/2), and optionally read name -> PS. Skips header / 'none'."""
    hp = {}
    ps = {} if with_ps else None
    n_none = 0
    with open(tsv) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            name, ht = p[0], p[1]
            if ht in ("H1", "1"):
                v = 1
            elif ht in ("H2", "2"):
                v = 2
            else:
                n_none += 1
                continue
            hp[name] = v
            if with_ps and len(p) >= 3 and p[2] not in ("none", ""):
                try:
                    ps[name] = int(p[2])
                except ValueError:
                    pass
    return hp, ps, n_none


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_bam")
    ap.add_argument("haplotags_tsv")
    ap.add_argument("out_bam")
    ap.add_argument("--threads", type=int, default=8, help="BAM (de)compression threads")
    ap.add_argument("--with-ps", action="store_true", help="also write the PS (phase set) tag")
    args = ap.parse_args()

    hp, ps, n_none = load_haplotags(args.haplotags_tsv, args.with_ps)
    print(f"[tag] reads with HP: {len(hp):,} | unphased (none): {n_none:,}", file=sys.stderr)

    n_in = n_tagged = n_stripped = 0
    with pysam.AlignmentFile(args.in_bam, "rb", threads=args.threads) as fin, \
            pysam.AlignmentFile(args.out_bam, "wb", template=fin, threads=args.threads) as fout:
        for read in fin:
            n_in += 1
            if read.has_tag("HP") or read.has_tag("PS"):
                n_stripped += 1
                read.set_tag("HP", None)
                read.set_tag("PS", None)
            v = hp.get(read.query_name)
            if v is not None:
                read.set_tag("HP", v, "i")
                if ps is not None and read.query_name in ps:
                    read.set_tag("PS", ps[read.query_name], "i")
                n_tagged += 1
            fout.write(read)

    print(f"[tag] records: {n_in:,} | tagged: {n_tagged:,} | pre-existing HP/PS removed: {n_stripped:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
