import ReadData as Data
import numpy as np
from collections import namedtuple

# Container for alignment result
AlignmentResult = namedtuple(
    'AlignmentResult',
    ['seq1', 'seq2', 'start1', 'start2',
     'end1', 'end2', 'n_gaps1', 'n_gaps2',
     'n_mismatches', 'score'])



class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self,  gap_start, gap_extend , blosum_mat, biosum_index):
        self.gap_open = gap_start
        self.gap_extend = gap_extend
        # 打分矩阵
        self.blosum_mat = blosum_mat
        self.blosum_index = biosum_index

    def get_two_char_score(self, a:chr,b:chr)->int:
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        return Data.get_Score_between_two_char(a,b,self.blosum_mat,self.blosum_index)

    def __str__(self):
        return "gap_start = %d; gap_extend = %d \n" % (
                self.gap_open, self.gap_extend
        )


class Sequence_Pair:
    def __init__(self,seq1:str = "" , seq2:str= ""):
        # 序列A，序列B
        if seq1 == "" and seq2 == "":
            pair = Data.Read_Two_Seq()
            self.sequence_A: str = pair[0]
            self.sequence_B: str = pair[1]
        elif seq1!="" and seq2 != "":
            self.sequence_A = seq1
            self.sequence_B = seq2
        # 统计结果
        self.gap_count = 0
        self.identity = 0
        self.similarity = 0
        self.final_length = 0



    def acc_identity(self):
        self.identity += 1

    def acc_gap_count(self):
        self.gap_count += 1

    def set_sequence_final_length(self, length):
        self.final_length = length

    def __str__(self):
        return "gap_count = %d; ideXntity = %d \n" % (self.gap_count, self.identity)



mat,index = Data.Read_BLO_Matrix()


def global_alignment(seqs = Sequence_Pair(), score = ScoreParam(-10,-0.5,mat,index),
                     method="global",max_hits=1):

    assert method == "global" or "local"

    seqi = seqs.sequence_A
    seqj = seqs.sequence_B

    NONE, LEFT, UP, DIAG = range(4)  # NONE is 0
    GAP_CHAR = '-'

    max_j = len(seqj)
    max_i = len(seqi)

    # 假定 seqi 比较长
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
    # 方向矩阵
    pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)  # NONE

    if method == 'global':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        F[0, 1:] = score.gap_open + score.gap_extend * \
            np.arange(0, max_j, dtype=np.float32)
        F[1:, 0] = score.gap_open + score.gap_extend * \
            np.arange(0, max_i, dtype=np.float32)

    for i in range(1, max_i + 1):
        ci = seqi[i-1: i]
        for j in range(1, max_j + 1):
            cj = seqj[j-1: j]
            # I
            I[i, j] = max(
                F[i, j - 1] + score.gap_open,
                I[i, j - 1] + score.gap_extend,
                J[i, j - 1] + score.gap_extend)
            # J
            J[i, j] = max(
                F[i - 1, j] + score.gap_open,
                J[i - 1, j] + score.gap_extend,
                I[i - 1, j] + score.gap_extend)
            # F
            diag_score = F[i - 1, j - 1] + score.get_two_char_score(ci,cj)
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
    ij_pairs = []
    #### !!!!!!!!!!!!!!!!!!!!!!
    ##   !!!!!!!!!!!!!!!!!!!!
    if method == 'local':
        # max anywhere
        maxv_indices = np.argwhere(F == F.max())[:max_hits]
        for index in maxv_indices:
            ij_pairs.append(index)
    elif method == "global":
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

        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
                ichar = seqi[i]
                jchar = seqj[j]
                if ichar != jchar:
                    n_mmatch += 1
                align_j.append(jchar)
                align_i.append(ichar)
            elif p == LEFT:
                j -= 1
                align_j.append(seqj[j])
                if not align_i or align_i[-1] != GAP_CHAR:
                    n_gaps_i += 1
                align_i.append(GAP_CHAR)
            elif p == UP:
                i -= 1
                align_i.append(seqi[i])
                if not align_j or align_j[-1] != GAP_CHAR:
                    n_gaps_j += 1
                align_j.append(GAP_CHAR)
            else:
                raise Exception('error!')
            p = pointer[i, j]
        align_i = ''.join(align_i)
        align_j = ''.join(align_j[::-1])

        aln = (AlignmentResult(align_i, align_j, i, j, end_i, end_j,
                               n_gaps_i, n_gaps_j, n_mmatch, score)
               if flip else
               AlignmentResult(align_j, align_i, j, i, end_j, end_i,
                               n_gaps_j, n_gaps_i, n_mmatch, score))

        results.append(aln)

    return results


seq1 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"+\
"KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLV"+\
"TLAAHLPAEFTP"+\
"AVHASLDKFLASVSTVLTSKYR"
seq2 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDL"+\
"STPDAVMGNPK"+\
"VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGN"+\
"VLVCVLAHHFG"+\
"KEFTPPVQAAYQKVVAGVANALAHKYH"

alg = global_alignment(seqs=Sequence_Pair(seq1,seq2))
print(alg)

