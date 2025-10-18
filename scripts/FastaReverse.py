#!/usr/bin/env python
# coding: utf-8

# In[33]:


import os
import textwrap

INPUT_DIR = "files"
OUT_DIR = os.path.join(INPUT_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

fasta_out_path = os.path.join(OUT_DIR, "revcomp.fasta")
summary_out_path = os.path.join(OUT_DIR, "summary.tsv")

def parse_fasta(path):
    """Генератор (header, sequence) для FASTA."""
    hdr = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # отдаем предыдущую запись
                if hdr is not None:
                    yield hdr, "".join(seq_chunks)
                hdr = line[1:].strip()  # сохраняем весь заголовок без '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # хвост
        if hdr is not None:
            yield hdr, "".join(seq_chunks)

def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

def wrap_fasta(seq, width=60):
    return "\n".join(textwrap.wrap(seq, width)) if width else seq

# Собираем список входных файлов
txt_files = [fn for fn in os.listdir(INPUT_DIR) if fn.endswith(".fasta")]

with open(fasta_out_path, "w") as fasta_out, open(summary_out_path, "w") as summary:
    summary.write("file\tseq_id\tlength\tinvalid_count\tinvalid_set\n")

    for fname in txt_files:
        in_path = os.path.join(INPUT_DIR, fname)
        for header, seq in parse_fasta(in_path):
            s = seq.upper()
            # считаем невалидные символы (не A/T/G/C/N)
            invalid = [ch for ch in s if ch not in "ATGCN"]
            rc = reverse_complement(s)

            # Пишем корректный FASTA
            fasta_out.write(f">{header}_RC\n")
            fasta_out.write(wrap_fasta(rc, 60) + "\n")

            # Пишем сводку
            inv_set = ",".join(sorted(set(invalid))) if invalid else "-"
            summary.write(f"{fname}\t{header}\t{len(s)}\t{len(invalid)}\t{inv_set}\n")

