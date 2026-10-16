"""Pin the headline findings so a change to the aligner cannot silently move them."""
from mthfr.analysis import run

res, _ = run()


def test_variant_residues_in_human():
    assert res["human"]["variants"]["222"]["human_residue"] == "A"
    assert res["human"]["variants"]["429"]["human_residue"] == "E"


def test_a222_conserved_in_all_vertebrates_and_maps_to_ecoli_a177():
    for acc in ["NP_034970.2", "NP_001315420.1", "XP_053257063.1", "NP_001096464.1", "NP_001268769.1"]:
        assert res["pairs"][acc]["residue_at"]["222"] == "A"
    assert res["pairs"]["NP_418376.1"]["residue_at"]["222"] == "A"
    assert res["deep_conservation"]["ecoli_residue_number_at_222"] == 177   # Guenther et al. 1999


def test_2023_fragment_cannot_see_either_variant():
    for acc in ["XP_028572270.1", "XP_028572269.1"]:
        assert res["pairs"][acc]["residue_at"] == {"222": "absent", "429": "absent"}
