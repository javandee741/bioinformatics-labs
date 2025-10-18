#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os

txt_files = [file for file in os.listdir("files/") if file.endswith('.fasta')]

print("file\tseq_id\tlength\tGC%")

for fname in txt_files:
    raw_dict = {}
    with open(os.path.join("files", fname)) as f:
        key = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]  # убираем '>' и берём id
                raw_dict[key] = ""
            else:
                raw_dict[key] += line.upper()

    for seq_id, seq in raw_dict.items():
        length = len(seq)
        if length == 0:
            gc_content = 0
        else:
            gc_content = 100 * (seq.count("C") + seq.count("G")) / length
        print(f"{fname}\t{seq_id}\t{length}\t{gc_content:.2f}")


# In[ ]:




