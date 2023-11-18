import pandas as pd
import sys
# motifs_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_positive/gene1/xstreme_DE/xstreme.txt"
# meme_motifs_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_positive/gene1/xstreme_DE/meme_out/meme.txt"

motifs_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_positive/"+sys.argv[1]+"/xstreme_DE/xstreme.txt"
meme_motifs_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_positive/"+sys.argv[1]+"/xstreme_DE/meme_out/meme.txt"

file_motifs = open(motifs_path, "r")
file_meme_motifs_path = open(meme_motifs_path, "r")

a = file_motifs.readlines()
b = file_meme_motifs_path.readlines()

temp_req_motif_list = []
for i in a:
    if "MEME" in i:
        if "MOTIF" in i:
            temp_req_motif_list.append(i.split())
# print(temp_req_motif_list)

req_motif_list = []
for i in range(0, len(temp_req_motif_list), 1):
    req_motif_list.append(temp_req_motif_list[i][1])
# print(req_motif_list)

l1 = []
for i in range(0, len(b), 1):
    for j in req_motif_list:
        if j and "BL  " in b[i]:
            for i in range(i, i+30, 1):
                if "//" in b[i]:
                    break
                else:
                    # print(b[i])
                    l1.append(b[i])
# print(l1)
indices = []
for i in range(0, len(l1), 1):
    if l1[i].startswith("BL"):
        indices.append(i)

l2 = []
for i in range(0, len(indices), 1):
    try:
        l2.append(l1[indices[i]:indices[i + 1]])
    except(IndexError):
        l2.append(l1[indices[i]:])
# print(l2)
l3 = []
for i in range(0, len(l2), 1):
    string1 = ""
    for j in range(0, len(l2[i]), 1):
        # string1 = ""
        string1 = string1 + l2[i][j]
    l3.append(string1)
# print(l3[6])
# print(l3)
# if "Seita" and"Sevir" and "Pahal" and "Pavir" and "Sobic" and "GRMZM" in l3[5]:
#     print("True")
# else:
#     print("False")

ans = []
for i in l3:
    # print(i)
    if "Seita" and"Sevir" and "Pahal" and "Pavir" and "Sobic" and "GRMZM" in i:
        # print("True")
        ans.append("True")
    else:
        # print("False")
        ans.append("False")
# print(len(ans))


name_motifs = []
for i in range(0, len(l2), 1):
    for j in range(0, len(l2[i]), 1):
        if l2[i][j].startswith("BL"):
            # print(l2[i][j])
            name_motifs.append(l2[i][j])

# print(len(name_motifs))

df = pd.DataFrame(list(zip(name_motifs, ans)), columns=["Motifs", "All_plants"])
print(df)
