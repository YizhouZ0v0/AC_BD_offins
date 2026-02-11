#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, gzip, sys
import pysam
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Filter reads by plasmid-aligned length/ratio using CIGAR from name-sorted BAM.")
    p.add_argument("--bam", required=True, help="Name-sorted BAM/SAM aligned to plasmid")
    p.add_argument("--r1", required=True, help="Input R1 FASTQ(.gz)")
    p.add_argument("--r2", help="Input R2 FASTQ(.gz); omit for single-end")
    p.add_argument("--out1", required=True, help="Output clean R1 FASTQ(.gz)")
    p.add_argument("--out2", help="Output clean R2 FASTQ(.gz); required for PE")
    p.add_argument("--min_match_len", type=int, default=50, help="Min total aligned length on query to plasmid to trigger removal [50]")
    p.add_argument("--min_match_frac", type=float, default=0.5, help="Min fraction of query aligned to plasmid to trigger removal [0.5]")
    p.add_argument("--drop_policy", choices=["either","both"], default="either",
                   help="either: if either mate exceeds threshold, drop both (recommended). both: require both mates exceed.")
    p.add_argument("--save_removed_prefix", help="If set, also write removed reads to {prefix}_R1/2.fastq.gz")
    return p.parse_args()

def merge_intervals(intervals):
    """intervals: list of (start, end) on query coords (0-based, end-exclusive).
       return total covered length after merging overlaps."""
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:  # overlap/adjacent
            cur_e = max(cur_e, e)
        else:
            total += (cur_e - cur_s)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s)
    return total

def collect_names_to_remove(bam_path, min_len, min_frac, drop_policy):
    """Scan name-sorted BAM; decide which read-names (R1/R2) to remove.
       Return two sets: names_R1_to_rm, names_R2_to_rm.
       For SE, only R1 set used.
    """
    bam = pysam.AlignmentFile(bam_path, "rb" if bam_path.endswith(".bam") else "r")
    names_R1_rm, names_R2_rm = set(), set()

    # group records by query_name (requires name-sorted)
    cur_name = None
    r1_intervals, r2_intervals = [], []
    r1_qlen, r2_qlen = None, None

    def decide_and_flush(name, r1_ints, r2_ints, r1_len, r2_len):
        if name is None:
            return
        # compute coverage and fraction per mate
        def decide(ints, qlen):
            if qlen is None or qlen == 0:
                return False
            cov = merge_intervals(ints)
            frac = cov / float(qlen)
            return (cov >= min_len) or (frac >= min_frac)

        r1_bad = decide(r1_ints, r1_len)
        r2_bad = decide(r2_ints, r2_len)

        if drop_policy == "either":
            bad_pair = r1_bad or r2_bad
        else:  # both
            # drop only if both mates exceed; for SE, this reduces to r1_bad
            bad_pair = r1_bad and (r2_bad if r2_len is not None else True)

        if bad_pair:
            names_R1_rm.add(name)
            if r2_len is not None:
                names_R2_rm.add(name)

    for aln in bam.fetch(until_eof=True):
        name = aln.query_name
        if name != cur_name:
            # flush previous
            decide_and_flush(cur_name, r1_intervals, r2_intervals, r1_qlen, r2_qlen)
            # reset
            cur_name = name
            r1_intervals, r2_intervals = [], []
            r1_qlen, r2_qlen = None, None

        if aln.is_unmapped:
            # unmapped contributes nothing; but we still want qlen for ratio
            if aln.is_read1:
                r1_qlen = max(r1_qlen or 0, aln.query_length or 0)
            elif aln.is_paired and aln.is_read2:
                r2_qlen = max(r2_qlen or 0, aln.query_length or 0)
            continue

        # query coords (0-based) covered by alignment (excludes soft-clipped bases)
        qs, qe = aln.query_alignment_start, aln.query_alignment_end  # pysam provides these
        if qs is None or qe is None:
            continue
        if aln.is_read1 or not aln.is_paired:
            r1_intervals.append((qs, qe))
            # prefer the max observed query_length across records
            ql = aln.query_length
            if ql is None:
                try:
                    ql = aln.infer_query_length(always=True)
                except Exception:
                    ql = None
            if ql is not None:
                r1_qlen = max(r1_qlen or 0, ql)
        else:  # read2
            r2_intervals.append((qs, qe))
            ql = aln.query_length
            if ql is None:
                try:
                    ql = aln.infer_query_length(always=True)
                except Exception:
                    ql = None
            if ql is not None:
                r2_qlen = max(r2_qlen or 0, ql)

    # flush last group
    decide_and_flush(cur_name, r1_intervals, r2_intervals, r1_qlen, r2_qlen)
    bam.close()
    return names_R1_rm, names_R2_rm

def fq_iter(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline()
            plus = fh.readline()
            q = fh.readline()
            yield h.rstrip(), s.rstrip(), plus.rstrip(), q.rstrip()

def fq_write(path, recs):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        for h,s,plus,q in recs:
            fh.write(f"{h}\n{s}\n{plus}\n{q}\n")

def header_to_name(header):
    """Normalize FASTQ header to read name (without /1 /2).
       Use first whitespace token; strip leading '@'; drop trailing /1,/2 if present."""
    tok = header[1:].split()[0]
    if tok.endswith("/1") or tok.endswith("/2"):
        tok = tok[:-2]
    return tok

def main():
    a = parse_args()
    pe = a.r2 is not None
    if pe and not a.out2:
        sys.exit("PE mode requires --out2")

    r1_rm, r2_rm = collect_names_to_remove(
        a.bam, a.min_match_len, a.min_match_frac, a.drop_policy
    )


    if pe and a.drop_policy == "either":
        union = r1_rm | r2_rm
        r1_rm = union
        r2_rm = union


    print(f"[INFO] Will remove {len(r1_rm)} reads from R1", file=sys.stderr)
    if pe:
        print(f"[INFO] Will remove {len(r2_rm)} reads from R2", file=sys.stderr)


    out_removed1 = out_removed2 = None
    if a.save_removed_prefix:
        out_removed1 = a.save_removed_prefix + "_R1.fastq.gz"
        if pe:
            out_removed2 = a.save_removed_prefix + "_R2.fastq.gz"

    def stream_filter(infq, outfq, rm_set, removed_out=None):
        def gen():
            for h,s,plus,q in fq_iter(infq):
                name = header_to_name(h)
                if name in rm_set:
                    if removed_out is not None:
                        removed_recs.append((h,s,plus,q))
                    continue
                yield (h,s,plus,q)

        if removed_out:
            removed_recs = []
            fq_write(outfq, gen())
            fq_write(removed_out, removed_recs)
        else:
            fq_write(outfq, gen())

    # R1
    stream_filter(a.r1, a.out1, r1_rm, out_removed1)
    # R2
    if pe:
        stream_filter(a.r2, a.out2, r2_rm, out_removed2)

if __name__ == "__main__":
    main()

