import math
import pandas as pd
# import os
from Bio import SeqIO

path_CDS = "/home/bgopikrishnan/Desktop/working_2/CDS/combined.fa"    # path to the directory that has all the CDS files
path_file = "/home/bgopikrishnan/Desktop/working_2/Putative_negative_ortholog_set/Putative_negative_ortholog_set.xlsx"

df1 = pd.read_excel(path_file)

# dropped columns that were not required for parsing.
df = df1.drop(["Orthogroup"], axis=1)
# print(df)
gene_dict = {}
for i in range(0, len(df), 1):
    gene_dict["list_%s" %str(i)] = []
# print(gene_dict)
for i in range(0, len(gene_dict), 1):
    gene_dict["list_%s" %str(i)] = list(df.iloc[i, :])
# print(gene_dict)
b = list(gene_dict.values())
# print(b[3])
# print(len(b))
seqs = []
for i in SeqIO.parse(path_CDS, "fasta"):
    seqs.append(i)
# print(len(seqs))
# print(b[0])
# trial_list = []
# for i in b[0]:
#     i = i.split(", ")
#     trial_list.extend(i)
# print(trial_list)

fl = []       # final list of genes
for i in range(0, len(b), 1):
    temp_list = []
    for j in b[i]:
        if type(j) == float:
            if math.isnan(j) == True:
                pass
        else:
            j = j.split(", ")
            temp_list.extend(j)
    fl.append(temp_list)

# print(len(fl))
# for i in fl:
#     print(len(i))
# genes_of_interest = []
# for seq in seqs:
    # for i in range(0, len(fl[4]), 1):
#         if (fl[4][i] in seq.description) == True:
#             genes_of_interest.append(seq)
# # print(len(genes_of_interest))

# SeqIO.write(genes_of_interest, "/home/bgopikrishnan/Desktop/working/muscle_input/try_2/input_4.fa", "fasta")


for i in range(0, len(fl), 1):
    genes_of_interest = []
    for j in fl[i]:
        # print(j)
        for k in seqs:
            # print(k.description)
            # print(type(k))
            if (j in k.description) == True:
                genes_of_interest.append(k)
    SeqIO.write(genes_of_interest, "/home/bgopikrishnan/Desktop/working_2/RAxML_input_putative_negative/gene_"+str(i+1)+".fa", "fasta")
# print(genes_of_interest)

