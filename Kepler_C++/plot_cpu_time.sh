#!/usr/bin/env bash
set -Eeuo pipefail

CSV="results/cpu_time/cpu_time_vs_h_all.csv"
OUT_DIR="results/plots/cpu_time"
OUT="${OUT_DIR}/cpu_time_kepler.pdf"
TITLE=""

PY="python3"
SCRIPT="./plotters/plot_cpu_time.py"

mkdir -p "${OUT_DIR}"

"${PY}" "${SCRIPT}" --csv "${CSV}" --out "${OUT}" --title "${TITLE}" --x h
