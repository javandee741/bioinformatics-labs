#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices


# ---------- utils ----------

def read_one_fasta(path: str):
    """Прочитать первую запись из FASTA: вернуть (id, seq_upper)."""
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        raise ValueError(f"No sequences in {path}")
    return recs[0].id, str(recs[0].seq).upper()


def load_blosum62_for_pairwise2() -> Dict[Tuple[str, str], int]:
    """
    Загрузить BLOSUM62 новым способом (Biopython 1.85+)
    и конвертировать в dict {(aa1, aa2): score} для pairwise2.globalds/localds.
    """
    mat = substitution_matrices.load("BLOSUM62")  # Matrix object (alphabet + numpy-матрица)
    aa = list(mat.alphabet)
    arr = np.array(mat)
    d = {}
    for i, a in enumerate(aa):
        for j, b in enumerate(aa):
            d[(a, b)] = int(arr[i, j])
    return d


def pick_matrix(name: str):
    """Вернуть матрицу замен для pairwise2 или None для 'simple' (+1/-1)."""
    if name.lower() == "blosum62":
        return load_blosum62_for_pairwise2()
    if name.lower() == "simple":
        return None
    raise ValueError("Unknown matrix (use 'blosum62' or 'simple')")


def do_align(seq1: str, seq2: str, method: str, matrix, gap_open: float, gap_extend: float):
    """
    Выполнить выравнивание:
      - если matrix is None: simple scoring (+1 / -1) через globalms/localms
      - иначе: dict-матрица через globalds/localds
    Вернёт лучший alignment из списка pairwise2 (a1, a2, score, start, end).
    """
    if matrix is None:
        # simple: match=+1, mismatch=-1
        if method == "global":
            alns = pairwise2.align.globalms(seq1, seq2, 1, -1, -gap_open, -gap_extend)
        else:
            alns = pairwise2.align.localms(seq1, seq2, 1, -1, -gap_open, -gap_extend)
    else:
        if method == "global":
            alns = pairwise2.align.globalds(seq1, seq2, matrix, -gap_open, -gap_extend)
        else:
            alns = pairwise2.align.localds(seq1, seq2, matrix, -gap_open, -gap_extend)
    if not alns:
        raise RuntimeError("Alignment failed: empty result")
    return alns[0]


def calc_metrics(a1: str, a2: str, score: float):
    """Посчитать length, matches, identity, gaps — по выровненным строкам с '-'."""
    assert len(a1) == len(a2)
    L = len(a1)
    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != "-" and y != "-")
    no_gaps_len = sum(1 for x, y in zip(a1, a2) if x != "-" and y != "-")
    identity = (matches / no_gaps_len * 100) if no_gaps_len else 0.0
    gaps1 = a1.count("-")
    gaps2 = a2.count("-")
    return {
        "length_aligned": L,
        "matches": matches,
        "no_gaps_len": no_gaps_len,
        "identity": identity,
        "gaps_seq1": gaps1,
        "gaps_seq2": gaps2,
        "score": score,
    }


def save_txt(out_txt: str, id1: str, id2: str, a1: str, a2: str, m: dict):
    os.makedirs(os.path.dirname(out_txt) or ".", exist_ok=True)
    with open(out_txt, "w", encoding="utf-8") as w:
        w.write("# Alignment report\n")
        w.write(f"# Seq1: {id1}\n# Seq2: {id2}\n\n")
        w.write(a1 + "\n" + a2 + "\n\n")
        for k, v in m.items():
            w.write(f"{k}\t{v}\n")


def save_tsv(out_tsv: str, id1: str, id2: str, m: dict):
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    header = [
        "seq1", "seq2", "score", "length_aligned",
        "no_gaps_len", "matches", "identity", "gaps_seq1", "gaps_seq2"
    ]
    row = [
        id1, id2, m["score"], m["length_aligned"],
        m["no_gaps_len"], m["matches"], f"{m['identity']:.2f}",
        m["gaps_seq1"], m["gaps_seq2"]
    ]
    new = not os.path.exists(out_tsv)
    with open(out_tsv, "a", encoding="utf-8") as w:
        if new:
            w.write("\t".join(header) + "\n")
        w.write("\t".join(map(str, row)) + "\n")


def save_png(out_png: str, identity: float):
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.figure(figsize=(4, 4))
    plt.bar(["Identity"], [identity], edgecolor="black")
    plt.ylim(0, 100)
    plt.ylabel("Percent (%)")
    plt.title("Alignment identity")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


# ---------- cli ----------

def main():
    ap = argparse.ArgumentParser(description="Align two sequences and report identity/score")
    ap.add_argument("-i1", required=True, help="FASTA 1")
    ap.add_argument("-i2", required=True, help="FASTA 2")
    ap.add_argument("-o", "--out", default="results/alignment.txt", help="output TXT")
    ap.add_argument("--tsv", default="results/alignment_summary.tsv", help="output TSV")
    ap.add_argument("--png", default="results/alignment_identity.png", help="output PNG")
    ap.add_argument("--method", choices=["global", "local"], default="global", help="alignment method")
    ap.add_argument("--matrix", choices=["blosum62", "simple"], default="blosum62", help="scoring scheme")
    ap.add_argument("--gap-open", type=float, default=10.0, help="gap opening penalty")
    ap.add_argument("--gap-extend", type=float, default=0.5, help="gap extension penalty")
    args = ap.parse_args()

    id1, s1 = read_one_fasta(args.i1)
    id2, s2 = read_one_fasta(args.i2)
    matrix = pick_matrix(args.matrix)
    a1, a2, score, start, end = do_align(s1, s2, args.method, matrix, args.gap_open, args.gap_extend)
    metrics = calc_metrics(a1, a2, score)

    save_txt(args.out, id1, id2, a1, a2, metrics)
    save_tsv(args.tsv, id1, id2, metrics)
    save_png(args.png, metrics["identity"])

    print(f"[OK] identity={metrics['identity']:.2f}% score={metrics['score']} -> "
          f"{os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()


# In[ ]:




