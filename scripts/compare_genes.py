#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys

def read_fasta(path: str) -> str:
    seq = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper().replace(" ", ""))
    return "".join(seq)


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)

def hamming(a: str, b: str) -> int:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i])

def main():
    ap = argparse.ArgumentParser(description="Compare two DNA sequences: lengths, GC%, Hamming")
    ap.add_argument("-i1", "--input1", required=True, help="FASTA 1 (DNA)")
    ap.add_argument("-i2", "--input2", required=True, help="FASTA 2 (DNA)")
    ap.add_argument("-o", "--out-tsv", required=True, help="Output TSV")
    args = ap.parse_args()

    s1 = read_fasta(args.input1)
    s2 = read_fasta(args.input2)

    L1, L2 = len(s1), len(s2)
    GC1, GC2 = gc_percent(s1), gc_percent(s2)
    hamm = hamming(s1, s2)

    with open(args.out_tsv, "w", encoding="utf-8") as w:
        w.write("Gene1\tGene2\tLength1\tLength2\tGC1\tGC2\tHamming\n")
        # имена возьмём из путей
        w.write(f"{args.input1}\t{args.input2}\t{L1}\t{L2}\t{GC1:.2f}\t{GC2:.2f}\t{hamm}\n")

    print(f"[OK] compare_genes → {args.out_tsv}")
    print(f"  Lengths: {L1} vs {L2}")
    print(f"  GC%%: {GC1:.2f} vs {GC2:.2f}")
    print(f"  Hamming (min length): {hamm}")

if __name__ == "__main__":
    sys.exit(main())
