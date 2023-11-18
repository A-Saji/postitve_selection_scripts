import sys
file_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_positive/"+sys.argv[1]+"/xstreme_combined/meme_out/meme.txt"
file = open(file_path, "r")
a = file.readlines()
# print(type(a))

l1 = []
for i in range(0, len(a), 1):
    # if "BL   " and "MOTIF" in a[i]:
    if "BL  " in a[i]:
        if "MOTIF" in a[i]:
            for j in range(i, i + 50, 1):
                if "//" in a[j]:
                    break
                else:
                    # print(a[j])
                    l1.append(a[j])

# print(l1)

indices = []
for i in range(0, len(l1), 1):
    if l1[i].startswith("BL"):
        indices.append(i)
# print(indices)

l2 = []
for i in range(0, len(indices), 1):
    try:
        l2.append(l1[indices[i]:indices[i + 1]])
    except(IndexError):
        l2.append(l1[indices[i]:])
# print(len(l2))

l3 = []
for i in range(0, len(l2), 1):
    string1 = ""
    for j in range(0, len(l2[i]), 1):
        # string1 = ""
        string1 = string1 + l2[i][j]
    l3.append(string1)


ans = []
for i in l3:
    if "Bradi" in i:
        continue
    if "LOC" in i:
        continue
    if "Brast" in i:
        continue
    if "Seita" in i:
        if "Sevir" in i:
            if "Pahal" in i:
                if "Sobic" in i:
                    if "GRMZM" in i:
                        ans.append(i)



ans_list = []
for i in ans:
    ans_list.append(i.split("\n"))

for i in ans_list:
    for j in i:
        if "BL" in j:
            print(j)