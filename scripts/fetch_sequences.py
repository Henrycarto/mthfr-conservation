"""Re-download data/mthfr_proteins.fasta from NCBI (needs network access).

    python scripts/fetch_sequences.py you@example.com
"""
import sys
from pathlib import Path

from Bio import Entrez, SeqIO

ACCESSIONS = {
    "NP_005948.3": "human", "NP_034970.2": "mouse", "NP_001315420.1": "chicken",
    "XP_053257063.1": "Aeolian wall lizard", "NP_001096464.1": "western clawed frog",
    "NP_001268769.1": "zebrafish",
    "XP_028572270.1": "common wall lizard, MTHFR-like (2023 assignment)",
    "XP_028572269.1": "common wall lizard, LOC114590064 isoform X1 (2023 assignment)",
    "NP_418376.1": "E. coli MetF",
}

if __name__ == "__main__":
    Entrez.email = sys.argv[1] if len(sys.argv) > 1 else "anonymous@example.com"
    out = Path(__file__).resolve().parent.parent / "data" / "fetched.fasta"
    with Entrez.efetch(db="protein", id=",".join(ACCESSIONS), rettype="fasta", retmode="text") as h:
        recs = list(SeqIO.parse(h, "fasta"))
    SeqIO.write(recs, out, "fasta")
    ref = {r.id: str(r.seq) for r in SeqIO.parse(out.parent / "mthfr_proteins.fasta", "fasta")}
    for r in recs:
        print(r.id, len(r.seq), "identical" if ref.get(r.id) == str(r.seq) else "DIFFERS")
