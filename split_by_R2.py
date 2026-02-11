#!/usr/bin/env python3
import sys, gzip

AC_SEQ = "TTCGGGAGAATTAGAGCTCCGGGTTACTGA"
BD_SEQ = "CTCGATGCAGACATGCTGCGACGATGATAT"

def open_gz(path, mode):
    return gzip.open(path, mode)

def write_record(w1, w2, r1, r2):
    w1.write(r1[0]); w1.write(r1[1]); w1.write(r1[2]); w1.write(r1[3])
    w2.write(r2[0]); w2.write(r2[1]); w2.write(r2[2]); w2.write(r2[3])

def main(r1_path, r2_path,prefix,only_record_target=True):
    # set only_record_target to speed up 
    outs = {}
    if only_record_target is True:
        for bin_name in ["AC","BD"]:
            outs[f"{bin_name}_R1"] = open_gz(f"tmp_{prefix}_{bin_name}_R1.fq.gz","wt")
            outs[f"{bin_name}_R2"] = open_gz(f"tmp_{prefix}_{bin_name}_R2.fq.gz","wt")
    else:
        for bin_name in ["AC","BD","other","conflict"]:
            outs[f"{bin_name}_R1"] = open_gz(f"tmp_{prefix}_{bin_name}_R1.fq.gz","wt")
            outs[f"{bin_name}_R2"] = open_gz(f"tmp_{prefix}_{bin_name}_R2.fq.gz","wt")

    total = ac = bd = other = conflict = 0

    with open_gz(r1_path,"rt") as r1, open_gz(r2_path,"rt") as r2:
        while True:
            r1_rec = [r1.readline() for _ in range(4)]
            r2_rec = [r2.readline() for _ in range(4)]

            if not r1_rec[0] or not r2_rec[0]:
                break
            total += 1


            seq2 = r2_rec[1].strip()

            has_ac = AC_SEQ in seq2
            has_bd = BD_SEQ in seq2

            if has_ac and not has_bd:
                write_record(outs["AC_R1"], outs["AC_R2"], r1_rec, r2_rec)
                ac += 1
            elif has_bd and not has_ac:
                write_record(outs["BD_R1"], outs["BD_R2"], r1_rec, r2_rec)
                bd += 1
            elif has_ac and has_bd:
                conflict += 1
                if only_record_target is True:
                    continue
                else:
                    write_record(outs["conflict_R1"], outs["conflict_R2"], r1_rec, r2_rec)
                
            else:
                other += 1
                if only_record_target is True:
                    continue
                else:
                    write_record(outs["other_R1"], outs["other_R2"], r1_rec, r2_rec)
                

    for f in outs.values():
        f.close()

    # summary
    print("Done.")
    print(f"Total:    {total}")
    print(f"AC:       {ac}")
    print(f"BD:       {bd}")
    print(f"Conflict: {conflict}")
    print(f"Other:    {other}")

if __name__ == "__main__":
    if len(sys.argv) !=4 :
        print(f"Usage: {sys.argv[0]} R1.fq.gz R2.fq.gz prefix", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2],sys.argv[3])

