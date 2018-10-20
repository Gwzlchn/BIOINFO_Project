from typing import Dict
import numpy as np

a1 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"+\
"KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLV"+\
"TLAAHLPAEFTP"+\
"AVHASLDKFLASVSTVLTSKYR"

def Read_Two_Seq(file_name: str = "Sequence.txt") -> list:
    Seq = ["", ""]
    f = open(file_name, "r")
    all_lines = f.readlines()
    # 读取序列
    Seq_cnt = -1
    for line in all_lines:
        if line[0] == '>':
            Seq_cnt += 1
            continue
        Seq[Seq_cnt] += line.strip('\n')
    f.close()
    assert Seq[0] == a1
    return Seq


def Read_BLO_Matrix(file_name: str = "BLOSUM62.txt") -> np.ndarray and dict:
    f = open(file_name, "r")
    all_lines = f.readlines()
    # Matrix 第一行是字母
    Matrix_Cnt = len(all_lines)
    Index = 0
    BLOSUM_Index_Dic: Dict[str, int] = {}
    for line in all_lines:
        if line[0] == '#':
            Matrix_Cnt -= 1
            continue
        elif line[0] == ' ':
            Matrix_Cnt -= 1
            BLOSUM_Matrix = np.zeros((Matrix_Cnt, Matrix_Cnt))
            continue
        elif line[0].isalpha() or line[0] == '*':
            # 映射关系建立
            BLOSUM_Index_Dic[line[0]] = Index
            # 矩阵建立
            one_line = line.strip('\n').split()
            one_line.remove(one_line[0])
            res = list(map(int, one_line))
            BLOSUM_Matrix[Index] = res
            Index += 1
    f.close()
    return BLOSUM_Matrix, BLOSUM_Index_Dic


def get_Score_between_two_char(a: chr, b: chr) -> int:
    mat, index = Read_BLO_Matrix("../BLOSUM62.txt")
    a_index = index[a.upper()]
    ret_list = mat[a_index]
    return ret_list[index[b.upper()]]

#print(get_Score_between_to_char('r','r'))
