# Looks for motifs which are exclusively present in C4 and not in C3 from the MEME.txt file in the meme directory of combined run.

import sys
file_path = "/home/bgopikrishnan/Desktop/new_working/XSTREME_input_putative_negative/"+sys.argv[1]+"/xstreme_combined/meme_out/meme.txt"
file = open(file_path, "r")
a = file.readlines()
# print(type(a))

l1 = []   # stores the block format into a list.
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

indices = []    # stores indices of the BL line from meme.txt file.
for i in range(0, len(l1), 1):
    if l1[i].startswith("BL"):
        indices.append(i)
# print(indices)

l2 = []     # extracts all the lines in between the two indices and makes a sub list in a main list.
for i in range(0, len(indices), 1):
    try:
        l2.append(l1[indices[i]:indices[i + 1]])
    except(IndexError):
        l2.append(l1[indices[i]:])       # Try except block catches the index error whichh will be thrown at the final iteration
# print(len(l2))

l3 = []      # sub lists are converted to strings and stored in list for easy string match search.
for i in range(0, len(l2), 1):
    string1 = ""
    for j in range(0, len(l2[i]), 1):
        # string1 = ""
        string1 = string1 + l2[i][j]
    l3.append(string1)


ans = []       # avoids considering motifs which is present in any C3 plants.
for i in l3:
    if "Bradi" in i:
        continue
    if "LOC" in i:
        continue
    if "Brast" in i:
        continue
    if "Seita" in i:        # considers motif if only present in all the C4 plants.
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
            print(j)   # prints the motifs which satisfy all these conditions.

    







