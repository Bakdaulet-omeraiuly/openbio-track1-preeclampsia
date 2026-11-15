"""POST HOC, descriptive: in the scramble arm of GSE80995, does 2% O2 (vs 20%) raise HIF2dep, hyp_all and FLT1?
Uses title-derived labels (see o1_amendment_title_labels.py). Usage: same args as o1_arnt.py.
"""
from pathlib import Path
import json
HERE = Path(__file__).resolve().parent
src = (HERE / "o1_arnt.py").read_text()
src = src.replace('info = pd.DataFrame({"sh": m["treatment"],',
                  'info = pd.DataFrame({"sh": title.str.extract(r"CTB_(Scramble|shARNT|Isolated)_")[0].map({"Scramble": "shScramble", "shARNT": "shARNT", "Isolated": "untreated"}),')
src = src[: src.index("res = {")]
ns = {"__file__": str(HERE / "o1_arnt.py"), "__name__": "o1"}
exec(compile(src, "o1_prefix", "exec"), ns)
z, cult, score, HYP, HIF2 = ns["z"], ns["cult"], ns["score"], ns["HYP"], ns["HIF2"]
out = {}
for name, s in [("hif2dep", score(HIF2)[0]), ("hyp_all", score(HYP)[0]), ("FLT1", z["FLT1"])]:
    d = {}
    for dn in sorted(cult.donor.unique()):
        sc = cult[(cult.donor == dn) & (cult.sh == "shScramble")]
        d[dn] = float(s[sc.index[sc.o2 == "hyp"]].iloc[0] - s[sc.index[sc.o2 == "norm"]].iloc[0])
    out[name] = {"mean_hyp_minus_norm_scramble": sum(d.values()) / len(d), "per_donor": d}
(HERE / "O1_posthoc_hypoxia_induction.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({k: round(v["mean_hyp_minus_norm_scramble"], 2) for k, v in out.items()}), {k: [round(x, 2) for x in v["per_donor"].values()] for k, v in out.items()})
