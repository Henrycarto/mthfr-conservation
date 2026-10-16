from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PURPLE, GREY, LIGHT, RED = "#53057F", "#9A9A9A", "#D9C7E6", "#C0392B"


def kmer_dotplot(a, b, k=3):
    """Fast word-match dot plot: a dot wherever a k-residue word in a occurs in b.
    Hash lookup, O(len(a) + len(b) + matches), instead of an O(n*m) aligner call per cell."""
    index = {}
    for j in range(len(b) - k + 1):
        index.setdefault(b[j:j + k], []).append(j)
    xs, ys = [], []
    for i in range(len(a) - k + 1):
        for j in index.get(a[i:i + k], ()):
            xs.append(i); ys.append(j)
    return np.array(xs), np.array(ys)


def make_figures(res, art, out):
    out = Path(out); out.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    seqs, cons = art["seqs"], np.array(art["cons"])
    h = seqs["NP_005948.3"]
    pairs = res["pairs"]

    # 1. conservation track along human MTHFR, coverage bars below
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 4.8), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1.3]))
    smooth = np.convolve(cons, np.ones(9) / 9, mode="same")
    x = np.arange(1, len(h) + 1)
    ax.fill_between(x, smooth, color=LIGHT, lw=0)
    ax.plot(x, smooth, color=PURPLE, lw=1)
    for p, lab in ((222, "A222V\n(C677T)"), (429, "E429A\n(A1298C)")):
        ax.axvline(p, color=RED, lw=1)
        ax.text(p + 6, 0.08, lab, color=RED, fontsize=8)
    ax.set_ylim(0, 1.05); ax.set_ylabel("identical in 5\northologs (9-aa mean)")
    ax.set_title("Human MTHFR: conservation across mouse, chicken, lizard, frog, zebrafish", loc="left")
    rows = [("NP_034970.2", "mouse"), ("XP_053257063.1", "lizard, P. raffonei"),
            ("NP_418376.1", "E. coli MetF"), ("XP_028572270.1", "P. muralis MTHFR-like (2023)")]
    for k, (acc, lab) in enumerate(rows):
        s, e = pairs[acc]["local_human_span"]
        col = RED if "2023" in lab else (GREY if "coli" in lab else PURPLE)
        ax2.barh(k, e - s + 1, left=s, color=col, height=0.6)
        ax2.text(len(h) + 8, k, lab, va="center", fontsize=8)
    for p in (222, 429):
        ax2.axvline(p, color=RED, lw=1)
    ax2.set_yticks([]); ax2.invert_yaxis(); ax2.set_xlim(0, len(h) + 175)
    ax2.set_xlabel("human MTHFR residue"); ax2.set_title("Region of human MTHFR each sequence aligns to (local alignment)", loc="left", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "conservation.png", dpi=160); plt.close(fig)

    # 2. dot plots: true lizard ortholog vs the 2023 'MTHFR-like' fragment
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
    for a, (acc, title) in zip(ax, [("XP_053257063.1", "vs P. raffonei MTHFR (656 aa)"),
                                    ("XP_028572270.1", "vs P. muralis MTHFR-like (219 aa), used in 2023")]):
        xs, ys = kmer_dotplot(h, seqs[acc], k=3)
        a.scatter(xs, ys, s=0.6, color=PURPLE)
        a.set_xlim(0, len(h)); a.set_ylim(len(seqs[acc]), 0)
        a.set_xlabel("human MTHFR residue"); a.set_ylabel("residue")
        a.set_title(title, loc="left", fontsize=9)
        a.axvline(222, color=RED, lw=0.8)
    fig.suptitle("Word-match dot plots (3-residue words). Red line: residue 222", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "dotplots.png", dpi=160); plt.close(fig)
