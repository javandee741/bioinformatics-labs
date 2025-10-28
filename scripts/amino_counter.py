#!/usr/bin/env python
# coding: utf-8

# In[9]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from collections import Counter
import matplotlib.pyplot as plt

STD_AA = list("ACDEFGHIKLMNPQRSTVWY")
SAFE_CHARS = re.compile(r"[^\w\-]+")

def _filter_jupyter_args():
    if any(arg == "-f" for arg in sys.argv):
        sys.argv = [sys.argv[0]]

def parse_fasta(path):
    header, chunks = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

def split_header(header):
    m = re.search(r"^(.*?)\s*\[(.+?)\]\s*$", header)
    if m:
        protein = m.group(1).strip() or "Unknown protein"
        organism = m.group(2).strip() or "Unknown organism"
    else:
        protein, organism = header, "Unknown organism"
    return protein, organism

def safe_name(s: str) -> str:
    return SAFE_CHARS.sub("_", s)[:120]

def count_aa(seq: str, collapse_unknown: bool = True):
    s = seq.upper()
    counts = Counter()
    for ch in s:
        if ch in STD_AA:
            counts[ch] += 1
        else:
            counts["X" if collapse_unknown else ch] += 1
    return counts, len(s)

def plot_bar(aa_counts, length: int, title: str, out_png: str):
    keys = [aa for aa in STD_AA if aa in aa_counts] + (["X"] if "X" in aa_counts else [])
    perc = [aa_counts[k] * 100.0 / max(1, length) for k in keys]
    plt.figure(figsize=(12, 6))
    bars = plt.bar(keys, perc, edgecolor="black")
    for b, p in zip(bars, perc):
        plt.text(b.get_x() + b.get_width() / 2.0, b.get_height() + 0.5,
                 f"{p:.1f}%", ha="center", va="bottom", fontsize=9)
    plt.title(title, fontsize=14)
    plt.xlabel("Amino acids", fontsize=12)
    plt.ylabel("Frequency (%)", fontsize=12)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    _filter_jupyter_args()

    ap = argparse.ArgumentParser(description="Amino-acid composition with TSV + bar plots")
    ap.add_argument("-i", "--input", default="files/protein.fasta", help="Входной FASTA-файл")
    ap.add_argument("-o", "--out-tsv", default="results/amino_stats.tsv", help="Выходной TSV-файл")
    ap.add_argument("--plots-dir", default="results", help="Папка для PNG-графиков")
    ap.add_argument("--min-length", type=int, default=20, help="Минимальная длина для графика")
    args = ap.parse_args()

    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for header, seq in parse_fasta(args.input):
        protein, organism = split_header(header)
        counts, length = count_aa(seq)

        total = max(1, length)
        for aa, cnt in sorted(counts.items()):
            rows.append([organism, protein, str(length), aa, str(cnt), f"{cnt * 100.0 / total:.2f}"])

        if length >= args.min_length:
            title = f"{protein} | {organism} (n={length})"
            png_name = f"amino_acid_freq_{safe_name(organism)}_{safe_name(protein)}.png"
            out_png = os.path.join(args.plots_dir, png_name)
            plot_bar(counts, length, title, out_png)
            print(f"[PNG] {os.path.abspath(out_png)}")

    with open(args.out_tsv, "w", encoding="utf-8") as w:
        w.write("organism\tprotein\tlength\taa\tcount\tpercent\n")
        for r in rows:
            w.write("\t".join(r) + "\n")

    print(f"[OK] TSV: {os.path.abspath(args.out_tsv)}")

if __name__ == "__main__":
    main()
