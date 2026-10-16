"""Pairwise protein alignment, written from scratch.

Needleman-Wunsch (global) and Smith-Waterman (local), both with affine gaps
(Gotoh's three-matrix recurrence). Gap convention matches EMBOSS and Biopython:
a gap of length L costs  open + (L - 1) * extend.

The substitution matrix (BLOSUM62) is data, loaded from Biopython; the dynamic
programming and traceback are implemented here and tested against Biopython's
PairwiseAligner in tests/.
"""
from dataclasses import dataclass

from Bio.Align import substitution_matrices

NEG = float("-inf")
BLOSUM62 = substitution_matrices.load("BLOSUM62")


@dataclass
class Alignment:
    a: str            # aligned first sequence, with '-'
    b: str            # aligned second sequence, with '-'
    score: float
    a_start: int      # 0-based start of the aligned region in the ungapped first sequence
    b_start: int

    @property
    def length(self):
        return len(self.a)

    @property
    def identities(self):
        return sum(x == y and x != "-" for x, y in zip(self.a, self.b))

    @property
    def identity(self):
        """EMBOSS convention: identical columns / alignment length (gaps included)."""
        return self.identities / self.length if self.length else 0.0

    def map_a_to_b(self):
        """For each 0-based position in sequence a covered by the alignment,
        the aligned residue in b ('-' if a gap)."""
        out, i = {}, self.a_start
        for x, y in zip(self.a, self.b):
            if x != "-":
                out[i] = y
                i += 1
        return out


def _align(s, t, local, gap_open, gap_extend, matrix):
    n, m = len(s), len(t)
    go, ge = gap_open, gap_extend   # negative numbers, e.g. -10, -0.5
    # M: ends in a match/mismatch, X: gap in t (consumes s), Y: gap in s (consumes t)
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0.0
    if local:
        for i in range(n + 1):
            M[i][0] = 0.0
        for j in range(m + 1):
            M[0][j] = 0.0
    else:
        for i in range(1, n + 1):
            X[i][0] = go + (i - 1) * ge
        for j in range(1, m + 1):
            Y[0][j] = go + (j - 1) * ge

    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        si = s[i - 1]
        Mi, Mp, Xi, Xp, Yi = M[i], M[i - 1], X[i], X[i - 1], Y[i]
        Yp = Y[i - 1]
        for j in range(1, m + 1):
            sub = matrix[si][t[j - 1]]
            diag = max(Mp[j - 1], Xp[j - 1], Yp[j - 1])
            mv = diag + sub
            if local and mv < 0:
                mv = 0.0
            Mi[j] = mv
            Xi[j] = max(Mp[j] + go, Xp[j] + ge, Yp[j] + go)
            Yi[j] = max(Mi[j - 1] + go, Yi[j - 1] + ge, Xi[j - 1] + go)
            if local and mv > best:
                best, bi, bj = mv, i, j

    # traceback
    if local:
        i, j, state, score = bi, bj, "M", best
    else:
        i, j = n, m
        finals = {"M": M[n][m], "X": X[n][m], "Y": Y[n][m]}
        state = max(finals, key=finals.get)
        score = finals[state]
    A, B = [], []
    while i > 0 or j > 0:
        if state == "M":
            if local and M[i][j] == 0:
                break
            sub = matrix[s[i - 1]][t[j - 1]]
            prev = M[i][j] - sub
            A.append(s[i - 1]); B.append(t[j - 1])
            i, j = i - 1, j - 1
            if i == 0 and j == 0:
                break
            state = next(k for k, v in (("M", M[i][j]), ("X", X[i][j]), ("Y", Y[i][j])) if v == prev)
        elif state == "X":
            v = X[i][j]
            A.append(s[i - 1]); B.append("-")
            i -= 1
            state = "X" if X[i][j] + ge == v else ("M" if M[i][j] + go == v else "Y")
        else:
            v = Y[i][j]
            A.append("-"); B.append(t[j - 1])
            j -= 1
            state = "Y" if Y[i][j] + ge == v else ("M" if M[i][j] + go == v else "X")
        if not local and (i == 0 or j == 0) and state == "M" and not (i == 0 and j == 0):
            state = "X" if j == 0 else "Y"
    return Alignment("".join(reversed(A)), "".join(reversed(B)), float(score), i, j)


def needleman_wunsch(s, t, gap_open=-10.0, gap_extend=-0.5, matrix=BLOSUM62):
    return _align(s, t, False, gap_open, gap_extend, matrix)


def smith_waterman(s, t, gap_open=-10.0, gap_extend=-0.5, matrix=BLOSUM62):
    return _align(s, t, True, gap_open, gap_extend, matrix)


def rescore(aln, gap_open=-10.0, gap_extend=-0.5, matrix=BLOSUM62):
    """Score an alignment from its columns; used to check traceback consistency."""
    total, prev_gap = 0.0, None
    for x, y in zip(aln.a, aln.b):
        if x == "-" or y == "-":
            kind = "a" if x == "-" else "b"
            total += gap_extend if prev_gap == kind else gap_open
            prev_gap = kind
        else:
            total += matrix[x][y]
            prev_gap = None
    return total
