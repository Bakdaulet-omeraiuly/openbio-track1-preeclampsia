"""Addendum O1 amendment (after the run): GSE80995 donor 630P normoxia arrays GSM2139930/GSM2139932 (corrected from a mistyped 'GSM2139927/GSM2139930' in the first version of this docstring) have
titles (Scramble / shARNT) that contradict their 'treatment' field. ARNT expression sides with the titles
(the treatment-labelled shARNT array had higher ARNT: delta_norm +1.82 vs -1.3 to -2.1 in other donors).
This reruns o1_arnt.py with shRNA taken from the sample TITLE. Hypoxia pairs are identical under both rules.
Usage: same as o1_arnt.py. Writes O1_results_title_labels.json.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
src = (HERE / "o1_arnt.py").read_text()
src = src.replace('info = pd.DataFrame({"sh": m["treatment"],',
                  'info = pd.DataFrame({"sh": title.str.extract(r"CTB_(Scramble|shARNT|Isolated)_")[0].map({"Scramble": "shScramble", "shARNT": "shARNT", "Isolated": "untreated"}),')
src = src.replace('"O1_results.json"', '"O1_results_title_labels.json"')
assert "O1_results_title_labels.json" in src and "CTB_(Scramble" in src
exec(compile(src, str(HERE / "o1_arnt.py"), "exec"), {"__file__": str(HERE / "o1_arnt.py"), "__name__": "__main__"})
