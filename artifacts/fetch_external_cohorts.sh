#!/usr/bin/env bash
# Download the public GEO files used by Addenda K, L, M and O into DIR (default: ./external_cohorts).
# Usage: bash artifacts/fetch_external_cohorts.sh [DIR]
# Then:
#   python artifacts/addendum-l/l1_hif2_program.py DIR/GSE216799_raw_counts.tsv.gz DIR/GSE216799_series_matrix.txt.gz
#   python artifacts/addendum-k/k_confirm_gse75010.py DIR/GSE75010_series_matrix.txt DIR/GPL6244.txt
#   python artifacts/addendum-l/l2_screen.py DIR/GSE75010_series_matrix.txt DIR/GPL6244.txt
#   python artifacts/addendum-l/l3_confirm.py DIR             # likewise l4_, l5_, l6_ (post hoc)
# Expected SHA-256 of the files we used are in artifacts/external_cohorts.sha256 (K, L) and
# artifacts/addendum-m/inputs.sha256 (M; run artifacts/addendum-m/m_confirm.py as described in its docstring).
set -euo pipefail
DIR="${1:-external_cohorts}"
mkdir -p "$DIR"
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series
get() { [ -s "$DIR/$2" ] || curl -fsSL --retry 3 -o "$DIR/$2" "$1"; echo "ok $2"; }
get "$GEO/GSE216nnn/GSE216799/suppl/GSE216799_raw_counts.tsv.gz" GSE216799_raw_counts.tsv.gz
get "$GEO/GSE216nnn/GSE216799/matrix/GSE216799_series_matrix.txt.gz" GSE216799_series_matrix.txt.gz
get "$GEO/GSE35nnn/GSE35574/matrix/GSE35574_series_matrix.txt.gz" GSE35574_series_matrix.txt.gz
get "$GEO/GSE24nnn/GSE24129/matrix/GSE24129_series_matrix.txt.gz" GSE24129_series_matrix.txt.gz
get "$GEO/GSE147nnn/GSE147776/matrix/GSE147776_series_matrix.txt.gz" GSE147776_series_matrix.txt.gz
get "$GEO/GSE75nnn/GSE75010/matrix/GSE75010_series_matrix.txt.gz" GSE75010_series_matrix.txt.gz
[ -s "$DIR/GSE75010_series_matrix.txt" ] || gunzip -k "$DIR/GSE75010_series_matrix.txt.gz"
# Addendum M (checksums in artifacts/addendum-m/inputs.sha256)
get "$GEO/GSE44nnn/GSE44711/matrix/GSE44711_series_matrix.txt.gz" GSE44711_series_matrix.txt.gz
get "$GEO/GSE74nnn/GSE74341/matrix/GSE74341_series_matrix.txt.gz" GSE74341_series_matrix.txt.gz
get "$GEO/GSE66nnn/GSE66273/matrix/GSE66273_series_matrix.txt.gz" GSE66273_series_matrix.txt.gz
# Addendum O1 (checksum in artifacts/addendum-o/inputs.sha256; uses GPL10558)
get "$GEO/GSE80nnn/GSE80995/matrix/GSE80995_series_matrix.txt.gz" GSE80995_series_matrix.txt.gz
for gpl in GPL6102 GPL6244 GPL20844 GPL10558 GPL16699 GPL4133; do
  get "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$gpl&targ=self&form=text&view=data" "$gpl.txt"
done
