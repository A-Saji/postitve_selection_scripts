from Bio import SeqIO
import sys
import pandas as pd

seqs = []
for i in SeqIO.parse("Dichanthelium_assembly.fna", "fasta"):
    seqs.append(i)

gene = sys.argv[1]
start_site = int(sys.argv[2]) - 1
# print(start_site)

for i in seqs:
    if (gene in i.description) == True:
        # print(i.seq)
        seq = i.seq
seq_string = str(seq)

substring = seq_string[0:start_site]
# print(substring)

upstream_seq = substring[-1: -1001 : -1]
# print(len(upstream_seq))

final_upstream = upstream_seq[-1: -1001: -1]
print(final_upstream.upper())
# print(len(final_upstream))