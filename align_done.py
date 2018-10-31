# -*- coding: utf-8 -*-

import sys


import numpy as np
from other_written.matrix import BLOSUM62
import ReadData as rd

__all__ = ["AlignmentResult", "aligner"]


class AlignmentResult:
    def __init__(self,seq1:str,seq2:str,seq_res:str,
                 start1:int, start2:int,end1:int,end2:int,
                 n_gap1:int, n_gap2:int,n_mismatches:int,
                 score:float):

        self.seq1 = seq1
        self.seq2 =seq2
        self.seq_res = seq_res
        self.start1 = start1
        self.start2 = start2
        self.end1 = end1
        self.end2 = end2
        self.n_gap1 = n_gap1
        self.n_gap2 = n_gap2
        self.n_mismatches = n_mismatches
        self.score = score

        #assert len(seq2) == len(seq1)
        #assert len(seq1) == len(self.seq_res)

    def __str__(self):
        final: str =""
        length = len(self.seq1)
        one_line_char = 50
        final += "SCORE:\t\t" + str(self.score) + "\n"
        final += "GAP1_COUNT:\t\t" + str(self.n_gap1) + "\n"
        final += "GAP2_COUNT:\t\t" + str(self.n_gap2) + "\n"
        final += "Start_1:\t\t" + str(self.start1) + "\n"
        final += "End_1:\t\t\t" + str(self.end1) + "\n"
        final += "Start_2:\t\t" + str(self.start2) + "\n"
        final += "End_2:\t\t\t" + str(self.end2) + "\n"
        #print(self.seq1)
        for i in range(int(length / one_line_char) + 1):
            i = i * one_line_char
            if i + one_line_char >= length:
                final +=  str(self.seq1[i: -1]) + "\n"
                final +=  str(self.seq_res[i:  -1])+ "\n"
                final +=  str(self.seq2[i: -1]) + "\n"
                break
            else:
                #print(seq1[i: i + one_line_char])
                final +=  str(self.seq1[i: i + one_line_char]) + "\n"
                final +=  str(self.seq_res[i:  i + one_line_char]) + "\n"
                final +=  str(self.seq2[i: i + one_line_char]) + "\n"
                final += "\n"

        return final



def aligner( seq_pair:list, matrix:np.ndarray , mat_index:dict,
             method:str='global',
            gap_open:float=-7, gap_extend:float=-7,
            gap_double=-7, max_hits=1):

    assert max_hits is None or max_hits > 0

    NONE, LEFT, UP, DIAG = range(4)  # NONE is 0
    GAP_CHAR = '-'
    seqi = seq_pair[0]
    seqj = seq_pair[1]

    max_j = len(seqj)
    max_i = len(seqi)
	
	#假定seq_i 为最长序列
    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i
    else:
        flip = 0

    F = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)
    I = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    I.fill(-np.inf)
    J = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    J.fill(-np.inf)
    pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)  # NONE

    if method == 'global':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        F[0, 1:] = gap_open + gap_extend * \
            np.arange(0, max_j, dtype=np.float32)
        F[1:, 0] = gap_open + gap_extend * \
            np.arange(0, max_i, dtype=np.float32)


    for i in range(1, max_i + 1):
        ci = seqi[i - 1:i]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1:j]
            # I
            I[i, j] = max(
                         F[i, j - 1] + gap_open,
                         I[i, j - 1] + gap_extend,
                         J[i, j - 1] + gap_double)
            # J
            J[i, j] = max(
                         F[i - 1, j] + gap_open,
                         J[i - 1, j] + gap_extend,
                         I[i - 1, j] + gap_double)
            # F
            #diag_score = F[i - 1, j - 1] + matrix[cj][ci]
            diag_score = F[i - 1, j - 1]  + rd.get_Score_between_two_char(cj,ci,matrix,mat_index)
            left_score = I[i, j]
            up_score = J[i, j]
            max_score = max(diag_score, up_score, left_score)

            F[i, j] = max(0, max_score) if method == 'local' else max_score

            if method == 'local':
                if F[i, j] == 0:
                    pass  # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i, j] = DIAG
                elif max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
            else:
                # global
                if max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
                else:
                    pointer[i, j] = DIAG

    # container for traceback coordinates
    ij_pairs = []

    if method == 'local':
        # max anywhere
        maxv_indices = np.argwhere(F == F.max())[:max_hits]
        for index in maxv_indices:
            ij_pairs.append(index)
    else:
        # method must be global at this point
        ij_pairs.append((i, j))

    results = []
    for i, j in ij_pairs:
        align_j = []
        align_i = []
        score = F[i, j]
        p = pointer[i, j]
        # mimic Python's coord system
        if method.startswith("global"):
            end_i, end_j = max_i, max_j
        else:
            end_i, end_j = i, j
        n_gaps_i, n_gaps_j, n_mmatch = 0, 0, 0
        seq_result = []
        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
                ichar = seqi[i]
                jchar = seqj[j]

                align_j.append(jchar)
                align_i.append(ichar)
                if ichar == jchar:
                    seq_result.append("|")
                else:
                    n_mmatch += 1
                    seq_result.append(" ")
            elif p == LEFT:
                seq_result.append("*")
                j -= 1
                align_j.append(seqj[j])
                if not align_i or align_i[-1] != GAP_CHAR:
                    n_gaps_i += 1
                align_i.append(GAP_CHAR)

            elif p == UP:
                seq_result.append("*")
                i -= 1
                align_i.append(seqi[i])
                if not align_j or align_j[-1] != GAP_CHAR:
                    n_gaps_j += 1
                align_j.append(GAP_CHAR)

            else:
                raise Exception('wtf!')
            p = pointer[i, j]

        align_i = ''.join(align_i[::-1])
        align_j = ''.join(align_j[::-1])
        seq_result_str = ''.join(seq_result[::-1])
        if flip:
            aln = AlignmentResult(align_i, align_j, seq_result_str , i, j, end_i, end_j,
                               n_gaps_i, n_gaps_j, n_mmatch, score)
        else:
            aln =  AlignmentResult(align_j, align_i, seq_result_str , j, i, end_j, end_i,
                               n_gaps_j, n_gaps_i, n_mmatch, score)

        results.append(aln)

    return results



# -----------------     test -------------------------------------------

mat, ind = rd.Read_BLO_Matrix(file_name="BLOSUM62.txt")
alg = aligner(rd.Read_Two_Seq(), mat, ind, gap_open=-10, gap_extend=-0.5,gap_double=-10,max_hits=2)
for i in alg:
    print(i)
alg=[]
alg = aligner(['ARAAV','ARAVVVARAV'],  mat, ind, method='local', gap_open=-10,gap_extend=-0.5,gap_double=-10,max_hits=2)
for i in alg:
    print(i)
alg = []
alg = aligner(rd.Read_Two_Seq(), mat, ind, method="local", gap_open=-10,gap_extend=-0.5,gap_double=-10,max_hits=5)
for i in alg:
    print(i)
    #print(len(i.seq1))
    print("----------------")