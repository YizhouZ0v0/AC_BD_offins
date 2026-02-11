#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import pysam
from collections import defaultdict
from typing import Optional, Tuple, List, Dict
import pandas as pd

def get_umi(aln: pysam.AlignedSegment, umi_tag: Optional[str], umi_regex: Optional[re.Pattern]) -> str:
    if umi_tag:
        try:
            return aln.get_tag(umi_tag)
        except KeyError:
            pass
    if umi_regex:
        m = umi_regex.search(aln.query_name or "")
        if m:
            return m.group(0)
    return "NOUMI"

def load_mapped_starts(path: str, umi_tag: Optional[str], umi_regex: Optional[re.Pattern]) -> Dict[str, List[Tuple[int, str, str]]]:
    data = defaultdict(list)  # chrom -> list of (start, umi, seq)
    with pysam.AlignmentFile(path, "r") as f:
        for aln in f.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            chrom = f.get_reference_name(aln.reference_id)
            start = aln.reference_start
            seq = aln.query_sequence or ""
            umi = get_umi(aln, umi_tag, umi_regex)
            data[chrom].append((start, umi, seq))
    return data

def cluster_hotspots(items: List[Tuple[int, str, str]], distance: int) -> List[Tuple[int, int, int, List[Tuple[int, str, str]]]]:
    if not items:
        return []
    items_sorted = sorted(items, key=lambda x: x[0])
    clusters, curr = [], [items_sorted[0]]
    for it in items_sorted[1:]:
        if it[0] - curr[-1][0] <= distance:
            curr.append(it)
        else:
            clusters.append(curr); curr = [it]
    clusters.append(curr)

    out = []
    for cluster in clusters:
        starts = [s for s,_,_ in cluster]
        start_min, start_max = min(starts), max(starts)
        center = (start_min + start_max) // 2
        out.append((start_min, start_max, center, cluster))
    return out

def summarize_clusters(per_chrom: Dict[str, List[Tuple[int, str, str]]], distance: int, sample_name: str) -> pd.DataFrame:
    rows = []
    for chrom, items in per_chrom.items():
        for start_min, start_max, center, cluster in cluster_hotspots(items, distance):
            n_reads = len(cluster)
            pair_set = set((umi, seq) for _, umi, seq in cluster)
            umi_set = set(umi for _, umi, _ in cluster)
            rows.append({
                "sample": sample_name, "chrom": chrom,
                "start_min": start_min, "start_max": start_max, "center": center,
                "n_reads": n_reads, "n_unique_pairs": len(pair_set), "n_unique_umis": len(umi_set),
            })
    cols = ["sample","chrom","start_min","start_max","center","n_reads","n_unique_pairs","n_unique_umis"]
    df = pd.DataFrame(rows, columns=cols)
    if len(df)==0:
        return df
    return df.sort_values(["chrom","center","start_min"]).reset_index(drop=True)

def pair_hotspots(dfA: pd.DataFrame, dfB: pd.DataFrame, pair_distance: int) -> pd.DataFrame:
    results = []
    for chrom in sorted(set(dfA["chrom"]).intersection(set(dfB["chrom"]))):
        a = dfA[dfA["chrom"]==chrom].sort_values("center").reset_index(drop=True)
        b = dfB[dfB["chrom"]==chrom].sort_values("center").reset_index(drop=True)
        i=j=0
        while i<len(a) and j<len(b):
            ca, cb = int(a.at[i,"center"]), int(b.at[j,"center"])
            diff = ca - cb
            if abs(diff) <= pair_distance:
                j_left=j
                while j_left-1>=0 and abs(ca-int(b.at[j_left-1,"center"]))<=pair_distance:
                    j_left-=1
                j_right=j
                while j_right+1<len(b) and abs(ca-int(b.at[j_right+1,"center"]))<=pair_distance:
                    j_right+=1
                for jj in range(j_left, j_right+1):
                    results.append({
                        "chrom": chrom,
                        "A_start_min": int(a.at[i,"start_min"]), "A_start_max": int(a.at[i,"start_max"]),
                        "A_center": ca, "A_n_reads": int(a.at[i,"n_reads"]),
                        "A_n_unique_pairs": int(a.at[i,"n_unique_pairs"]), "A_n_unique_umis": int(a.at[i,"n_unique_umis"]),
                        "B_start_min": int(b.at[jj,"start_min"]), "B_start_max": int(b.at[jj,"start_max"]),
                        "B_center": int(b.at[jj,"center"]), "B_n_reads": int(b.at[jj,"n_reads"]),
                        "B_n_unique_pairs": int(b.at[jj,"n_unique_pairs"]), "B_n_unique_umis": int(b.at[jj,"n_unique_umis"]),
                        "center_distance": abs(ca - int(b.at[jj,"center"])),
                    })
                if ca <= cb: i+=1
                else: j+=1
            elif diff < 0:
                i+=1
            else:
                j+=1
    cols = [
        "chrom",
        "A_start_min","A_start_max","A_center","A_n_reads","A_n_unique_pairs","A_n_unique_umis",
        "B_start_min","B_start_max","B_center","B_n_reads","B_n_unique_pairs","B_n_unique_umis",
        "center_distance",
    ]
    df = pd.DataFrame(results, columns=cols)
    if len(df)==0:
        return df
    return df.sort_values(["chrom","A_center","B_center","center_distance"]).reset_index(drop=True)

def main():
    ap = argparse.ArgumentParser(description="Cluster mapped-start hotspots with UMI+sequence dedup and joint A/B pairing.")
    ap.add_argument("--a", required=True, help="Sample A SAM/BAM path")
    ap.add_argument("--b", required=True, help="Sample B SAM/BAM path")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--cluster-distance", type=int, default=1000, help="Distance (bp) to merge starts into a hotspot (default: 1000)")
    ap.add_argument("--pair-distance", type=int, default=1000, help="Max center distance (bp) to pair A/B hotspots (default: 1000)")
    ap.add_argument("--umi-tag", default=None, help="SAM tag for UMI (e.g., RX or UB). If set, overrides --umi-regex")
    ap.add_argument("--umi-regex", default="(?<=:UMI_)[A-Za-z0-9_-]+", help="Regex to extract UMI from read name (default matches ':UMI_XXXX').")
    args = ap.parse_args()

    umi_regex = re.compile(args.umi_regex) if args.umi_regex else None

    a_map = load_mapped_starts(args.a, args.umi_tag, umi_regex)
    b_map = load_mapped_starts(args.b, args.umi_tag, umi_regex)

    dfA = summarize_clusters(a_map, args.cluster_distance, sample_name="A")
    dfB = summarize_clusters(b_map, args.cluster_distance, sample_name="B")

    out_prefix = args.out_prefix
    a_out = f"{out_prefix}.A.hotspots.tsv"
    b_out = f"{out_prefix}.B.hotspots.tsv"
    dfA.to_csv(a_out, sep="\t", index=False)
    dfB.to_csv(b_out, sep="\t", index=False)

    pair_df = pair_hotspots(dfA, dfB, args.pair_distance)
    pair_out = f"{out_prefix}.AB.pairs.tsv"
    pair_df.to_csv(pair_out, sep="\t", index=False)

    summ_df = pd.DataFrame([{
        "A_n_hotspots": len(dfA),
        "B_n_hotspots": len(dfB),
        "AB_n_pairs": len(pair_df),
        "cluster_distance_bp": args.cluster_distance,
        "pair_distance_bp": args.pair_distance,
    }])
    summ_df.to_csv(f"{out_prefix}.summary.tsv", sep="\t", index=False)

    print("Done.")
    print(f"A hotspots: {a_out}")
    print(f"B hotspots: {b_out}")
    print(f"AB pairs  : {pair_out}")
    print(f"Summary   : {out_prefix}.summary.tsv")

if __name__ == "__main__":
    main()

