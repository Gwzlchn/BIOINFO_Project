import numpy as np

Seq = ["",""]
f = open("Sequence.txt","r")  
all_lines = f.readlines()
# 读取序列
Seq_cnt = -1
for line in all_lines:
    if line[0] == '>':
        Seq_cnt+=1
        continue
    Seq[Seq_cnt] += line.strip('\n')
   
f.close()

f = open("BLOSUM62.txt","r")
all_lines = f.readlines() 
# Matrix 第一行是字母
Matrix_Cnt = len(all_lines)
BLOSUM_INDEX = []
for line in all_lines:
    if line[0] == '#' :
        Matrix_Cnt -= 1
        continue
    if line[0] == ' ':
        Matrix_Cnt -= 1
        BLOSUM_Matrix = np.zeros((Matrix_Cnt, Matrix_Cnt))
        
        continue
    
    
print(Matrix_Cnt)
f.close()