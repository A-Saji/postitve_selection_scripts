# Name: - B. Gopikrishnan
# Date: - 9th March 2023

# Code to get the PACIDs into a txt file which can be submitted to Phytozome to get upstream seqs
import pandas as pd
from Bio import SeqIO

CDS_file = "/home/bgopikrishnan/Desktop/working/CDS/Zmays.fa"   # Have to provide the CDS file.
df = pd.read_excel("test.xlsx")  # A temp excel file from which the column can be read. (copy paste one column)
# print(df)
a = df["Zmays_284_Ensembl-18_2010-01-MaizeSequence.protein_primaryTranscriptOnly"]    # provide the column name

gene_list = []
for i in a:
    if str(i) == "nan":
        continue
    else:
        gene_list.append(i)

trial_list = []
for i in gene_list:
    trial_list.append(i.split(", "))
# print(trial_list)

gl = []
for i in range(0, len(trial_list), 1):
    for j in trial_list[i]:
        gl.append(j)
# print(gl)
seqs = []
for i in SeqIO.parse(CDS_file, "fasta"):
    seqs.append(i)


description_list = []
for i in gl:
    for j in range(0, len(seqs), 1):
        if i in seqs[j].description:
            description_list.append(seqs[j].description)
        else:
            pass

description_list_splitted = []
for i in description_list:
    description_list_splitted.append(i.split(" "))


pacid_list = []
for i in range(0, len(description_list_splitted), 1):
    pacid_list.append(description_list_splitted[i][1])

temp_pacid_list = []
for i in pacid_list:
    temp_pacid_list.append(i.split("="))


final_pacid_list = []
for i in range(0, len(temp_pacid_list), 1):
    final_pacid_list.append(temp_pacid_list[i][1])
# print(len(final_pacid_list))

file = open("Zmays.txt", "w")   # Writes the output of all PACIDs into a txt file.
for i in final_pacid_list:
    file.write(i+"\n")
file.close()



