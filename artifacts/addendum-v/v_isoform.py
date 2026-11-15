"""Addendum V: EPAS1-null vs HIF1A-null effects on FLT1, the HIF2dep program and Hallmark hypoxia in iPSC-derived
primitive syncytium at 3% O2 (rules in plan.md).

  python v_isoform.py fetch DIR
  python v_isoform.py run DIR
"""

import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S = pd.read_csv(HERE / "samples.tsv", sep="\t")
URL = {"GSE288288": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288288/suppl/GSE288288_study1_genesCounts.csv.gz",
       "GSE302997": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302997/suppl/GSE302997_study8_genesCounts.csv.gz"}
TARGET = {"GSE288288": "EPAS1", "GSE302997": "HIF1A"}
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    for g, u in URL.items():
        out = d / u.rsplit("/", 1)[1]
        if not out.exists():
            urllib.request.urlretrieve(u, out)
        print("ok", out.name)


COLS = {
    "GSE288288": {"null_hyp": ["KO_A03__Hypoxia", "KO_B01__Hypoxia", "KO_B02__Hypoxia", "CE_E06__Hypoxia", "CE_G07__Hypoxia",
                               "CE_H06__Hypoxia", "PTC_A09__Hypoxia", "PTC_D10__Hypoxia", "PTC_G09__Hypoxia"],
                  "null_norm": ["KO_A03__Normoxia", "KO_B01__Normoxia", "KO_B02__Normoxia", "CE_E06__Normoxia", "CE_G07__Normoxia",
                                "CE_H06__Normoxia", "PTC_A09__Normoxia", "PTC_D10__Normoxia", "PTC_G09__Normoxia"],
                  "wt_hyp": ["WT_1__Hypoxia", "WT_2__Hypoxia", "WT_3__Hypoxia"],
                  "wt_norm": ["WT_1__Normoxia", "WT_2__Normoxia", "WT_3__Normoxia"]},
    "GSE302997": {"null_hyp": ["HIF1a_PTC_A04_3_O2", "HIF1a_PTC_A05_3_O2", "HIF1a_PTC_A06_3_O2"],
                  "null_norm": ["HIF1a_PTC_A04_20_O2", "HIF1a_PTC_A05_20_O2", "HIF1a_PTC_A06_20_O2"],
                  "wt_hyp": ["KOLF2_2J_3_O2_REP1", "KOLF2_2J_3_O2_REP2", "KOLF2_2J_3_O2_REP3"],
                  "wt_norm": ["KOLF2_2J_20_O2_REP1", "KOLF2_2J_20_O2_REP2", "KOLF2_2J_20_O2_REP3"]},
}
STRATEGY = {"KO": "Full coding length", "CE": "Critical exon", "PTC": "Termination codon"}


def load(d, study):
    c = pd.read_csv(d / URL[study].rsplit("/", 1)[1], sep="\t", index_col=0)
    c = c.select_dtypes(include=[np.number])
    rows = []
    for group, prefixes in COLS[study].items():
        for pre in prefixes:
            hits = [col for col in c.columns if col.startswith(pre + "_")]
            assert len(hits) == 1, (pre, hits)
            strategy = STRATEGY.get(pre.split("_")[0], "Termination codon") if group.startswith("null") else "Not applicable"
            rows.append({"col": hits[0], "gene": TARGET[study] if group.startswith("null") else "Not applicable",
                         "strategy": strategy, "condition": "hypoxia" if group.endswith("hyp") else "normoxia"})
    meta = pd.DataFrame(rows)
    assert meta.col.is_unique
    c = c[meta.col.tolist()]
    m = pd.read_csv(ROOT / "data_external/ensembl_gene_id_to_name_GRCh37.87.tsv", sep="\t", header=None, names=["id", "sym"])
    m = m.drop_duplicates("id").set_index("id")["sym"]
    ids = c.index.astype(str).str.split(".").str[0]
    sym = pd.Series(ids, index=c.index).map(m)
    c = c[sym.notna().values].groupby(sym.dropna().values).sum()
    assert {"FLT1", "EPAS1", "HIF1A"} <= set(c.index)
    return c, meta


def welch_less(a, b):
    """one-sided p for mean(a) < mean(b), and the difference a - b"""
    t = stats.ttest_ind(a, b, equal_var=False, alternative="less")
    return float(np.mean(a) - np.mean(b)), float(t.pvalue)


def welch_greater(a, b):
    t = stats.ttest_ind(a, b, equal_var=False, alternative="greater")
    return float(np.mean(a) - np.mean(b)), float(t.pvalue)


def analyse(d, study):
    c, meta = load(d, study)
    cpm = c / c.sum() * 1e6
    keep = (cpm >= 1).sum(axis=1) >= 3
    lg = np.log2(cpm[keep] + 1)
    z = lg.sub(lg.mean(1), axis=0).div(lg.std(1).replace(0, np.nan), axis=0)
    score = lambda genes: z.loc[z.index.intersection(sorted(genes))].mean()
    tab = pd.DataFrame({"program": score(HIF2), "hyp_all": score(HYP),
                        "FLT1": lg.loc["FLT1"] if "FLT1" in lg.index else np.nan,
                        "target": lg.loc[TARGET[study]] if TARGET[study] in lg.index else np.nan}, index=lg.columns)
    tab = tab.join(meta.set_index("col")[["gene", "strategy", "condition"]])
    tab["null"] = tab.gene != "Not applicable"
    wt_h = tab[~tab.null & (tab.condition == "hypoxia")]; wt_n = tab[~tab.null & (tab.condition == "normoxia")]
    ko_h = tab[tab.null & (tab.condition == "hypoxia")]; ko_n = tab[tab.null & (tab.condition == "normoxia")]
    out = {"n": {"wt_hyp": len(wt_h), "wt_norm": len(wt_n), "null_hyp": len(ko_h), "null_norm": len(ko_n)},
           "n_genes_program": int(len(z.index.intersection(sorted(HIF2)))), "n_genes_hyp": int(len(z.index.intersection(sorted(HYP))))}
    out["gate1_wt_hypoxia_raises_hyp_all"] = dict(zip(("diff", "p"), welch_greater(wt_h.hyp_all, wt_n.hyp_all)))
    out["gate1_passed"] = bool(out["gate1_wt_hypoxia_raises_hyp_all"]["diff"] > 0 and out["gate1_wt_hypoxia_raises_hyp_all"]["p"] < 0.05)
    out["gate2_target_mrna_null_vs_wt_hyp"] = dict(zip(("diff", "p"), welch_less(ko_h.target, wt_h.target)))
    out["wt_hypoxia_effect_on_FLT1"] = dict(zip(("diff", "p_two_sided"), (float(wt_h.FLT1.mean() - wt_n.FLT1.mean()),
                                            float(stats.ttest_ind(wt_h.FLT1, wt_n.FLT1, equal_var=False).pvalue))))
    for v in ("FLT1", "program", "hyp_all"):
        out[f"null_vs_wt_at_3pct_{v}"] = dict(zip(("diff", "p_one_sided_lower"), welch_less(ko_h[v], wt_h[v])))
        # interaction: (null_h - null_n) - (wt_h - wt_n)
        out[f"interaction_{v}"] = float((ko_h[v].mean() - ko_n[v].mean()) - (wt_h[v].mean() - wt_n[v].mean()))
    if study == "GSE288288":
        fl = ko_h[ko_h.strategy == "Full coding length"]
        out["sens_a_full_deletion_only"] = {v: dict(zip(("diff", "p_one_sided_lower"), welch_less(fl[v], wt_h[v])))
                                            for v in ("FLT1", "program", "hyp_all")}
    out["per_sample"] = tab.reset_index().rename(columns={"index": "column"}).to_dict(orient="records")
    return out


def call(e2, e1):
    """isoform reading from plan.md, given (diff, p) for EPAS1-null and HIF1A-null."""
    s2, s1 = e2[1] < 0.05 and e2[0] < 0, e1[1] < 0.05 and e1[0] < 0
    if s2 and s1:
        return "both"
    if s2 and abs(e1[0]) < 0.5 * abs(e2[0]):
        return "HIF-2-preferring"
    if s1 and abs(e2[0]) < 0.5 * abs(e1[0]):
        return "HIF-1-preferring"
    return "neither / unresolved"


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    res = {g: analyse(d, g) for g in URL}
    e, h = res["GSE288288"], res["GSE302997"]
    reading = {}
    for v in ("FLT1", "program"):
        if e["gate1_passed"] and h["gate1_passed"]:
            k = f"null_vs_wt_at_3pct_{v}"
            reading[v] = call((e[k]["diff"], e[k]["p_one_sided_lower"]), (h[k]["diff"], h[k]["p_one_sided_lower"]))
        else:
            reading[v] = "not testable (gate 1 failed)"
    res["isoform_reading"] = reading
    (HERE / "V_results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    for g, r in ((k, res[k]) for k in URL):
        print("==", g, TARGET[g], r["n"], "genes prog/hyp", r["n_genes_program"], r["n_genes_hyp"])
        print("  gate1", r["gate1_wt_hypoxia_raises_hyp_all"], "passed", r["gate1_passed"])
        print("  gate2 target mRNA", r["gate2_target_mrna_null_vs_wt_hyp"])
        print("  WT hypoxia effect on FLT1", r["wt_hypoxia_effect_on_FLT1"])
        for v in ("FLT1", "program", "hyp_all"):
            print(f"  null vs WT @3% {v}", r[f"null_vs_wt_at_3pct_{v}"], "interaction %.3f" % r[f"interaction_{v}"])
        if "sens_a_full_deletion_only" in r:
            print("  sens full deletion", r["sens_a_full_deletion_only"])
    print("ISOFORM READING", reading)


if __name__ == "__main__":
    main()
