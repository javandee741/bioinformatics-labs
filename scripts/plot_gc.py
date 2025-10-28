#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, csv
import matplotlib.pyplot as plt
import os

def main():
    ap = argparse.ArgumentParser(description="Plot GC%% comparison from TSV")
    ap.add_argument("-i", "--input", required=True, help="compare_genes.tsv")
    ap.add_argument("-o", "--out-png", required=True, help="output PNG")
    args = ap.parse_args()

    with open(args.input, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        row = next(r)  # одна строка
        labels = ["seq1", "seq2"]
        values = [float(row["GC1"]), float(row["GC2"])]



    plt.figure(figsize=(6,4))
    bars = plt.bar(labels, values, edgecolor="black")
    for b, v in zip(bars, values):
        plt.text(b.get_x() + b.get_width()/2, v + 0.5, f"{v:.1f}%", ha="center", va="bottom")
    plt.ylabel("GC%")
    plt.title("GC comparison")
    plt.ylim(0, 100)
    plt.tight_layout()
    os.makedirs(os.path.dirname(args.out_png) or ".", exist_ok=True)
    plt.savefig(args.out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[OK] plot → {args.out_png}")

if __name__ == "__main__":
    main()

