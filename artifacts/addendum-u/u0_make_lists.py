"""Addendum U, step 0: split the HIF2dep program into ARNT-bound / unbound genes using Addendum T's rule (±10 kb of an NM_ TSS).
Uses ChIP and annotation only (no tissue data). Usage: python u0_make_lists.py T_DATA_DIR
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "artifacts/addendum-t"))
import t_chip as T  # noqa: E402

peaks, ref = T.load(Path(sys.argv[1]))
bound = T.bound_genes(peaks, ref, 10_000)
expr = pd.read_csv(ROOT / "artifacts/addendum-l/L1_GSE216799_PT2385_vs_DMSO.csv", index_col=0)["baseMean"]
universe = set(expr[expr >= 10].index) & set(ref.gene)
genes = [g for g in (ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split() if g in universe]
here = Path(__file__).resolve().parent
(here / "HIF2dep_ARNT_bound.txt").write_text("\n".join(sorted(g for g in genes if g in bound)) + "\n")
(here / "HIF2dep_ARNT_unbound.txt").write_text("\n".join(sorted(g for g in genes if g not in bound)) + "\n")
