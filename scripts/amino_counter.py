#!/usr/bin/env python
# coding: utf-8

# In[9]:


import os
import textwrap
import re
import matplotlib.pyplot as plt
import numpy as np



raw_dict = {} 
with open("../data_samples/protein.fasta") as f:
    current_key = None
    current_value = None
    
    for line in f:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith(">"):
            # Используем re.search вместо re.findall для более надежного извлечения
            pattern1 = r'>(.*?) \['  # название белка до квадратной скобки
            pattern2 = r'\[(.*?)\]'   # содержимое в квадратных скобках
            
            match1 = re.search(pattern1, line)
            match2 = re.search(pattern2, line)
            
            if match1 and match2:
                key = match2.group(1)   # содержимое в квадратных скобках (организм)
                value = match1.group(1) # название белка
                
                # Инициализируем вложенный словарь, если ключа еще нет
                if key not in raw_dict:
                    raw_dict[key] = {}
                
                # Сохраняем текущие ключи для последующих строк с последовательностью
                current_key = key
                current_value = value
                
                # Инициализируем пустую строку для последовательности белка
                if current_value not in raw_dict[current_key]:
                    raw_dict[current_key][current_value] = ""
                    
        else:
            # Добавляем последовательность к текущему белку
            if current_key is not None and current_value is not None:
                raw_dict[current_key][current_value] += line

print("Структура данных:")
for organism, proteins in raw_dict.items():
    print(f"Организм: {organism}")
    for protein_name, sequence in proteins.items():
        amino_dict = {}
        for i in sequence:
            key_to_check = i
            if key_to_check in amino_dict:
                amino_dict[key_to_check] += 1
            else:
                amino_dict[key_to_check] = 1
        print(f"\tБелок: {protein_name}")
        print(f"\tПоследовательность: {sequence[:50]}...")  # показываем первые 50 символов
        print (f"\tАминокислоты:")
        for key, value in amino_dict.items():
            print(f"\t\t{key} {100*value/len(sequence):.2f}")
        print(f"\tДлина: {len(sequence)}")
        print()


# Создаем графики для каждого белка
for organism, proteins in raw_dict.items():
    for protein_name, sequence in proteins.items():
        # Считаем частоты аминокислот
        amino_dict = {}
        for i in sequence:
            if i in amino_dict:
                amino_dict[i] += 1
            else:
                amino_dict[i] = 1
        
        # Сортируем аминокислоты по частоте
        sorted_amino = sorted(amino_dict.items(), key=lambda x: x[1], reverse=True)
        amino_acids = [item[0] for item in sorted_amino]
        frequencies = [item[1] for item in sorted_amino]
        percentages = [100 * freq / len(sequence) for freq in frequencies]
        
        # Создаем график
        plt.figure(figsize=(12, 6))
        bars = plt.bar(amino_acids, percentages, color='skyblue', edgecolor='black')
        
        # Добавляем значения на столбцы
        for bar, percentage in zip(bars, percentages):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{percentage:.1f}%', ha='center', va='bottom', fontsize=9)
        
        plt.title(f'Частоты аминокислот: {protein_name}\nОрганизм: {organism}', fontsize=14)
        plt.xlabel('Аминокислоты', fontsize=12)
        plt.ylabel('Частота (%)', fontsize=12)
        plt.xticks(rotation=45)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        # Сохраняем график
        safe_protein_name = re.sub(r'[^\w\-_]', '_', protein_name)
        safe_organism = re.sub(r'[^\w\-_]', '_', organism)
        plt.savefig(f'amino_acid_freq_{safe_organism}_{safe_protein_name}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График сохранен для: {protein_name} ({organism})")

print("Все графики созданы!")


    






# In[22]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from collections import Counter

import matplotlib.pyplot as plt

# 20 стандартных аминокислот
STD_AA = list("ACDEFGHIKLMNPQRSTVWY")
SAFE_CHARS = re.compile(r"[^\w\-]+")

def _filter_jupyter_args():
    """Убираем служебные аргументы Jupyter (-f <path>), чтобы argparse не падал."""
    # Пример: ['amino_counter.py', '-f', 'C:\\...\\kernel-xxxx.json']
    if any(arg == "-f" for arg in sys.argv):
        sys.argv = [sys.argv[0]]

def parse_fasta(path):
    """Генератор (header, seq) из FASTA; игнорируем пустые строки."""
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
    """
    Пробуем извлечь Protein и Organism из 'Protein [Organism]'.
    Если нет квадратных скобок — возвращаем (header, 'Unknown organism').
    """
    m = re.search(r"^(.*?)\s*\[(.+?)\]\s*$", header)
    if m:
        protein = m.group(1).strip() or "Unknown protein"
        organism = m.group(2).strip() or "Unknown organism"
    else:
        protein, organism = header, "Unknown organism"
    return protein, organism

def safe_name(s: str) -> str:
    """Безопасное имя файла (латинские буквы/цифры/подчёрки/дефис)."""
    return SAFE_CHARS.sub("_", s)[:120]

def count_aa(seq: str, collapse_unknown: bool = True):
    """Подсчёт частот AA; всё вне STD_AA -> 'X' (если collapse_unknown=True)."""
    s = seq.upper()
    counts = Counter()
    for ch in s:
        if ch in STD_AA:
            counts[ch] += 1
        else:
            counts["X" if collapse_unknown else ch] += 1
    return counts, len(s)

def plot_bar(aa_counts, length: int, title: str, out_png: str):
    """Сохраняем бар-чарт частот AA в PNG."""
    # Фиксируем порядок столбцов: STD_AA (+ 'X', если есть)
    keys = [aa for aa in STD_AA if aa in aa_counts] + (["X"] if "X" in aa_counts else [])
    perc = [aa_counts[k] * 100.0 / max(1, length) for k in keys]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(keys, perc, edgecolor="black")
    for b, p in zip(bars, perc):
        plt.text(
            b.get_x() + b.get_width() / 2.0,
            b.get_height() + 0.5,
            f"{p:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    plt.title(title, fontsize=14)
    plt.xlabel("Amino acids", fontsize=12)
    plt.ylabel("Frequency (%)", fontsize=12)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    _filter_jupyter_args()  # чтобы работало из Jupyter без ошибок argparse

    ap = argparse.ArgumentParser(
        description="Amino-acid composition with TSV + per-protein bar plots"
    )
    ap.add_argument(
        "-i", "--input",
        default="../data_samples/protein.fasta",
        help="Входной FASTA (по умолчанию: ../data_samples/protein.fasta)",
    )
    ap.add_argument(
        "-o", "--out-tsv",
        default="../data_samples/results_amino/amino_stats.tsv",
        help="Выходной TSV (по умолчанию: ../data_samples/results_amino/amino_stats.tsv)",
    )
    ap.add_argument(
        "--plots-dir",
        default="results",
        help="Каталог для PNG-графиков (по умолчанию: results)",
    )
    ap.add_argument(
        "--min-length",
        type=int,
        default=20,
        help="Минимальная длина последовательности для построения графика (по умолчанию: 20)",
    )
    ap.add_argument(
        "--keep-unknown",
        action="store_true",
        help="Не сворачивать нестандартные аминокислоты в 'X'",
    )
    args = ap.parse_args()

    if args.plots_dir == "results" and args.out_tsv:
    # если пользователь не менял plots-dir,
    # по умолчанию складываем PNG туда же, куда и TSV
        tsv_dir = os.path.dirname(args.out_tsv) or "."
        args.plots_dir = tsv_dir

    # Готовим каталог для TSV
    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    n_records = 0
    for header, seq in parse_fasta(args.input):
        n_records += 1
        protein, organism = split_header(header)
        counts, length = count_aa(seq, collapse_unknown=not args.keep_unknown)

        # Накапливаем TSV-строки
        total = max(1, length)
        for aa, cnt in sorted(counts.items()):
            rows.append([
                organism,
                protein,
                str(length),
                aa,
                str(cnt),
                f"{cnt * 100.0 / total:.2f}",
            ])

        # Сохраняем график (если длина достаточно большая)
        if length >= args.min_length:
            title = f"{protein} | {organism} (n={length})"
            png_name = f"amino_acid_freq_{safe_name(organism)}_{safe_name(protein)}.png"
            out_png = os.path.join(args.plots_dir, png_name)
            plot_bar(counts, length, title, out_png)

    # Пишем TSV
    with open(args.out_tsv, "w", encoding="utf-8") as w:
        w.write("organism\tprotein\tlength\taa\tcount\tpercent\n")
        for r in rows:
            w.write("\t".join(r) + "\n")

    print(f"[OK] Прочитано записей: {n_records}")
    print(f"[OK] TSV: {args.out_tsv}")
    print(f"[OK] PNG-схемы (если были): {os.path.abspath(args.plots_dir)}")

if __name__ == "__main__":
    main()


# In[ ]:




