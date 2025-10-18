#!/usr/bin/env python
# coding: utf-8

# In[ ]:


raw_dict 
with open("files/example.fasta") as f:
    key = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            key = line
            raw_dict[key] = ""
        else:
            raw_dict[key] += line.upper()  # на всякий случай приводим к верхнему регистру

for key, seq in raw_dict.items():
    stringHeader = key.split()[0]
    length = len(seq)
    gc_content = 100 * (seq.count("C") + seq.count("G")) / length
    invalid = [b for b in seq if b not in "ATGCN"]

    print(stringHeader)
    print(f"Length: {length}")
    print(f"GC%: {gc_content:.2f}")
    print(f"Invalid symbols: {len(invalid)} ({set(invalid) if invalid else '-'})\n")


# In[ ]:




