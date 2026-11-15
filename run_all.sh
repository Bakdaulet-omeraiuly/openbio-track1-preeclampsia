#!/usr/bin/env bash
# One command to reproduce everything: notebooks 01-12, every addendum script (downloads public GEO files into
# ./data_cache, checksummed per addendum), and the summary figures.
# Usage: bash run_all.sh            (about 1-2 hours and ~2 GB of downloads; notebooks alone take ~10 min)
#        bash run_all.sh notebooks  (notebooks + summary figures only; no downloads beyond the datathon package)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python}
(cd notebooks && for nb in 01_*.ipynb 02_*.ipynb 03_*.ipynb 05_*.ipynb 06_*.ipynb 07_*.ipynb 08_*.ipynb 09_*.ipynb 10_*.ipynb 11_*.ipynb 12_*.ipynb; do
  jupyter nbconvert --to notebook --execute --inplace "$nb"; done)
if [ "${1:-all}" != "notebooks" ]; then
  C=data_cache; mkdir -p $C
  bash artifacts/fetch_external_cohorts.sh $C/ext
  $PY artifacts/addendum-k/k_confirm_gse75010.py $C/ext/GSE75010_series_matrix.txt $C/ext/GPL6244.txt
  $PY artifacts/addendum-l/l1_hif2_program.py $C/ext/GSE216799_raw_counts.tsv.gz $C/ext/GSE216799_series_matrix.txt.gz
  $PY artifacts/addendum-l/l3_confirm.py $C/ext
  for s in l4_calibration_posthoc l5_stb_adjustment_posthoc l6_deconvolution_posthoc l8_diagnosis_sensitivity_posthoc; do $PY artifacts/addendum-l/$s.py $C/ext; done
  $PY artifacts/addendum-m/m_confirm.py $C/ext
  $PY artifacts/addendum-n/n_explore.py $C/ext/GSE75010_series_matrix.txt $C/ext/GPL6244.txt
  $PY artifacts/addendum-o/o1_arnt.py $C/ext/GSE80995_series_matrix.txt.gz $C/ext/GPL10558.txt
  $PY artifacts/addendum-o2/o2_within_celltype.py
  G=https://ftp.ncbi.nlm.nih.gov/geo/series
  curl -fsSL -o $C/GSE186257_series_matrix.txt.gz $G/GSE186nnn/GSE186257/matrix/GSE186257_series_matrix.txt.gz
  curl -fsSL -o $C/GSE186257_norm.txt.gz $G/GSE186nnn/GSE186257/suppl/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz
  curl -fsSL -o $C/GSE324038_counts.txt.gz $G/GSE324nnn/GSE324038/suppl/GSE324038_Project_Knm32.rsem_genes_expected_count.txt.gz
  $PY artifacts/addendum-p/p_confirm.py $C/GSE186257_series_matrix.txt.gz $C/GSE186257_norm.txt.gz
  $PY artifacts/addendum-q/q_hif1_vs_hif2.py $C/GSE324038_counts.txt.gz
  for a in "r/r_within_donor r" "s/s_spatial s" "t/t_chip t" "v/v_isoform v" "w/w_mvm w" "x/x_infarct x"; do
    set -- $a; $PY artifacts/addendum-$1.py fetch $C/$2; $PY artifacts/addendum-$1.py run $C/$2; done
  $PY artifacts/addendum-r/r3_purity_posthoc.py $C/r
  $PY artifacts/addendum-u/u0_make_lists.py $C/t && $PY artifacts/addendum-u/u_direct_targets.py $C/s
  $PY artifacts/addendum-w/w2_calibration_posthoc.py $C/w
  $PY artifacts/addendum-y/y_maternal.py $C/ext $C/w
fi
$PY artifacts/summary-figures/make_figures.py
echo "done"
