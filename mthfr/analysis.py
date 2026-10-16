"""Is the residue hit by MTHFR C677T conserved across vertebrates?

    python -m mthfr.analysis      # writes results/results.json and figures/

Human MTHFR (NP_005948.3, 656 aa) is aligned to its orthologs in five
vertebrates, to the two Podarcis muralis "MTHFR-like" proteins used in the 2023
course assignment, and to E. coli MetF. Every alignment uses the from-scratch
aligners in mthfr/align.py (BLOSUM62, gap open -10, extend -0.5) and is
cross-checked against Biopython's score at run time.
"""
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

from .align import BLOSUM62, needleman_wunsch, smith_waterman

ROOT = Path(__file__).resolve().parent.parent
HUMAN = "NP_005948.3"
VERTEBRATES = ["NP_034970.2", "NP_001315420.1", "XP_053257063.1", "NP_001096464.1", "NP_001268769.1"]
ASSIGNMENT_2023 = ["XP_028572270.1", "XP_028572269.1"]
ECOLI = "NP_418376.1"

VARIANTS = {   # human 1-based residue -> (variant, dbSNP, cDNA change)
    222: ("A222V", "rs1801133", "c.665C>T, historically C677T"),
    429: ("E429A", "rs1801131", "c.1286A>C, historically A1298C"),
}


def load():
    recs = {r.id: r for r in SeqIO.parse(ROOT / "data" / "mthfr_proteins.fasta", "fasta")}
    return {k: str(v.seq) for k, v in recs.items()}, {k: v.description.split(" ", 1)[1] for k, v in recs.items()}


def check_against_biopython(aln, s, t, mode):
    a = PairwiseAligner()
    a.substitution_matrix = BLOSUM62
    a.open_gap_score, a.extend_gap_score, a.mode = -10.0, -0.5, mode
    ref = a.score(s, t)
    assert abs(ref - aln.score) < 1e-6, (mode, ref, aln.score)


def _col(aln, a_pos):
    """Alignment column index of 0-based position a_pos in the first sequence."""
    i = aln.a_start
    for col, x in enumerate(aln.a):
        if x != "-":
            if i == a_pos:
                return col
            i += 1
    raise KeyError(a_pos)


def window(seq_map, center, half=8):
    return "".join(seq_map.get(i, " ") for i in range(center - half - 1, center + half))


def run():
    seqs, names = load()
    h = seqs[HUMAN]
    res = {"human": {"accession": HUMAN, "length": len(h),
                     "variants": {str(p): dict(variant=v[0], dbsnp=v[1], cdna=v[2], human_residue=h[p - 1])
                                  for p, v in VARIANTS.items()}},
           "pairs": {}}
    alignments = {}
    for acc in VERTEBRATES + ASSIGNMENT_2023 + [ECOLI]:
        s = seqs[acc]
        g = needleman_wunsch(h, s)
        loc = smith_waterman(h, s)
        check_against_biopython(g, h, s, "global")
        check_against_biopython(loc, h, s, "local")
        use = loc if acc in ASSIGNMENT_2023 + [ECOLI] else g
        m = use.map_a_to_b()
        covered = [i for i, c in m.items()]
        human_span = (min(covered) + 1, max(covered) + 1)
        res["pairs"][acc] = dict(
            name=names[acc], length=len(s),
            global_identity=round(g.identity, 3), global_score=g.score,
            local_identity=round(loc.identity, 3), local_score=loc.score,
            local_human_span=[loc.a_start + 1, loc.a_start + len(loc.a.replace("-", ""))],
            human_coverage=round(len(loc.a.replace("-", "")) / len(h), 3),
            residue_at={str(p): m.get(p - 1, "absent") for p in VARIANTS},
            window_222=window(m, 222),
        )
        alignments[acc] = (g, loc)
    res["human"]["window_222"] = window({i: c for i, c in enumerate(h)}, 222)

    # per-position conservation across the five vertebrate orthologs (human coordinates)
    maps = [alignments[a][0].map_a_to_b() for a in VERTEBRATES]
    cons = [sum(mp.get(i) == h[i] for mp in maps) / len(maps) for i in range(len(h))]
    res["conservation"] = dict(
        mean_identity_to_human=round(sum(cons) / len(cons), 3),
        fully_conserved_positions=int(sum(c == 1.0 for c in cons)),
        at_222=cons[221], at_429=cons[428],
        # how unusual is full conservation? share of positions conserved in all five
        share_fully_conserved=round(sum(c == 1.0 for c in cons) / len(cons), 3))

    # stricter test: inside the region E. coli MetF aligns to, which positions are
    # identical in all five vertebrates AND in the bacterial enzyme?
    em = alignments[ECOLI][1].map_a_to_b()
    region = sorted(em)
    deep = [i for i in region if cons[i] == 1.0 and em[i] == h[i]]
    res["deep_conservation"] = dict(
        ecoli_aligned_human_span=[region[0] + 1, region[-1] + 1],
        positions_in_region=len(region),
        identical_vertebrates_and_ecoli=len(deep),
        share=round(len(deep) / len(region), 3),
        residue_222_is_one_of_them=(221 in deep),
        ecoli_residue_number_at_222=alignments[ECOLI][1].b_start + 1 + sum(
            1 for x, y in zip(alignments[ECOLI][1].a[:_col(alignments[ECOLI][1], 221)], alignments[ECOLI][1].b[:_col(alignments[ECOLI][1], 221)]) if y != "-"))
    return res, dict(seqs=seqs, alignments=alignments, cons=cons)


def main():
    res, art = run()
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "results.json").write_text(json.dumps(res, indent=2))
    from .figures import make_figures
    make_figures(res, art, ROOT / "figures")
    print(f"human 222 = {res['human']['variants']['222']['human_residue']}, 429 = {res['human']['variants']['429']['human_residue']}")
    print(f"{'accession':16s} {'name':55s} {'len':>4s} {'gID':>5s} {'lID':>5s} {'cov':>5s} {'span':>10s}  222 429  window_222")
    print(f"{'':16s} {'':55s} {'':4s} {'':5s} {'':5s} {'':5s} {'':10s}           {res['human']['window_222']}  (human)")
    for acc, p in res["pairs"].items():
        print(f"{acc:16s} {p['name'][:55]:55s} {p['length']:4d} {p['global_identity']:5.2f} {p['local_identity']:5.2f} "
              f"{p['human_coverage']:5.2f} {str(p['local_human_span']):>10s}  {p['residue_at']['222']:>3s} {p['residue_at']['429']:>3s}  {p['window_222']}")
    print(res["conservation"])
    print(res["deep_conservation"])


if __name__ == "__main__":
    main()
