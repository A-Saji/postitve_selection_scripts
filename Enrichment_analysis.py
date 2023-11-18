import pandas as pd
import sys
import os
# file = open("meme.txt", "r")
file = open("/home/bgopikrishnan/Desktop/working_2/XSTREME_input_negative/"+sys.argv[1]+"/xstreme_combined/meme_out/meme.txt")
a = file.readlines()
# print(a)
l1 = []
for i in range(0, len(a), 1):
    # if "BL   " and "MOTIF" in a[i]:
    if "BL  " in a[i]:
        if "MOTIF" in a[i]:
            for j in range(i, i + 50, 1):
                if "//" in a[j]:
                    l1.append('//')
                    break
                else:
                    # print(a[j])
                    l1.append(a[j])

C4_set = ["ELECO", "Misin", "Oropetium", "Pavir", "Seita", "Urofu", "GRMZM"]
C3_set = ["Bradi", "Brast", "Chala", "lcl", "HORVU", "LOC"]

# print(l1)
size = len(l1)
idx_list = [idx + 1 for idx, val in
            enumerate(l1) if val == "//"]
res = [l1[i: j] for i, j in
        zip([0] + idx_list, idx_list +
        ([size] if idx_list[-1] != size else []))]

motif_name_list = []
for i in range(0, len(res), 1):
    motif_name_list.append(res[i][0])
# print(motif_name_list)
motif_name_list_stripped = []

for i in motif_name_list:
    motif_name_list_stripped.append(i.rstrip())
# print(motif_name_list_stripped)

l3 = []
for i in res:
    str1 = ""
    l3.append(str1.join(i))
# print(l3[4])

# C4_counter = 0
# C3_counter = 0

C4_vals = []
for i in l3:
    C4_counter = 0
    for j in C4_set:
        if j in i:
            C4_counter += 1
    C4_vals.append(C4_counter)
# print(C4_vals)

C3_vals = []
for i in l3:
    C3_counter = 0
    for j in C3_set:
        if j in i:
            C3_counter += 1
    C3_vals.append(C3_counter)

# print(C3_vals)

C4_total_number = []
C3_total_number = []

for i in range(0, len(l3), 1):
    C4_total_number.append(7)
    C3_total_number.append(6)

C4_percentage = []
for i in C4_vals:
    C4_percentage.append((i/7)*100)

C3_percentage = []
for i in C3_vals:
    C3_percentage.append((i/6)*100)

delta = []
for i in range(0, len(C4_percentage), 1):
    delta.append(C4_percentage[i] - C3_percentage[i])

df = pd.DataFrame(list(zip(motif_name_list_stripped, C4_vals, C4_total_number, C4_percentage, C3_vals, C3_total_number, C3_percentage, delta)), columns=["Motif", "C4 occ", "Total C4", "C4_perc", "C3 occ", "Total C3", "C3_perc", "Delta"])
# print(df)

enriched_list = []
for i in range(0, len(delta), 1):
    if delta[i] >= 50:
        # print(motif_name_list[i])
        enriched_list.append(motif_name_list_stripped[i])
print(enriched_list)

# df.to_excel(sys_argv[1]+"_enrichment_analysis.xlsx")
# df.to_excel("enrichment_analysis_"+sys.argv[1]+".xlsx")
# os.system("mv enrichment_analysis_"+sys.argv[1]+".xlsx /home/bgopikrishnan/Desktop/working_2/XSTREME_input_negative/"+sys.argv[1])

# with open("/home/bgopikrishnan/Desktop/working_2/XSTREME_input_negative/"+sys.argv[1]+"/"+sys.argv[1]+"_enriched.txt", "w") as f:
#     for i in enriched_list:
#         f.write(i)
