import random

import pytest
from Bio.Align import PairwiseAligner

from mthfr.align import BLOSUM62, needleman_wunsch, rescore, smith_waterman


def bio(mode):
    a = PairwiseAligner()
    a.substitution_matrix = BLOSUM62
    a.open_gap_score, a.extend_gap_score, a.mode = -10.0, -0.5, mode
    return a


def rand_protein(rng, n):
    return "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


def mutate(rng, s, rate=0.25):
    out = []
    for c in s:
        r = rng.random()
        if r < rate * 0.6:
            out.append(rng.choice("ACDEFGHIKLMNPQRSTVWY"))
        elif r < rate * 0.8:
            continue                                    # deletion
        elif r < rate:
            out.append(c + rand_protein(rng, rng.randint(1, 4)))  # insertion
        else:
            out.append(c)
    return "".join(out)


@pytest.mark.parametrize("seed", range(25))
def test_scores_match_biopython(seed):
    rng = random.Random(seed)
    s = rand_protein(rng, rng.randint(20, 120))
    t = mutate(rng, s) if seed % 3 else rand_protein(rng, rng.randint(20, 120))
    g, l = needleman_wunsch(s, t), smith_waterman(s, t)
    assert g.score == pytest.approx(bio("global").score(s, t))
    assert l.score == pytest.approx(bio("local").score(s, t))
    # traceback reproduces the score and the input sequences
    assert rescore(g) == pytest.approx(g.score)
    assert rescore(l) == pytest.approx(l.score)
    assert g.a.replace("-", "") == s and g.b.replace("-", "") == t
    assert s[l.a_start:].startswith(l.a.replace("-", ""))
    assert t[l.b_start:].startswith(l.b.replace("-", ""))


def test_identical_sequences():
    s = "MVNEARGNSSLNPCLEG"
    g = needleman_wunsch(s, s)
    assert g.identity == 1.0 and g.a == s and g.b == s


def test_position_mapping():
    g = needleman_wunsch("ACDEFG", "ACEFG")       # D deleted in b
    m = g.map_a_to_b()
    assert m[0] == "A" and m[2] == "-" and m[5] == "G"
