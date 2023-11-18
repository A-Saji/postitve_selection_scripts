# raxmlHPC-PTHREADS-SSE3 -f a -m PROTGAMMAAUTO -p 10000 -x 10000 -s input_1_aligned.aln -# 200 -T 3 -n input_1_N_tree.txt
# muscle -align input_1_for_muscle.fa -output input_1_aligned.aln

import os

path = "/home/bgopikrishnan/Desktop/working/muscle_alignment"
file_names = os.listdir(path)
file_names.sort()
# print(file_names)

for i in range(0, len(file_names), 1):
    os.system("raxmlHPC-PTHREADS-SSE3 -f a -m PROTGAMMAAUTO -p 10000 -x 10000 -s /home/bgopikrishnan/Desktop/working/muscle_alignment/Alignment_{}/aln_{}.aln -# 200 -T 3 -n txt".format(i, i))