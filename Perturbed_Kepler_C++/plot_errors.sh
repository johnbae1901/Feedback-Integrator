#!/usr/bin/env bash
set -Eeuo pipefail

# ---- Inputs produced by your C++ run ----
ERRORS_CSV="results/error_data/errors_vs_h.csv"
MAXDL_CSV="results/error_data/maxdL_vs_h.csv"
MAXDA_CSV="results/error_data/maxdA_vs_h.csv"

# ---- Output (empty => show interactively) ----
OUT="results/plots/error/error_kepler.pdf"   # e.g., .png/.pdf; set "" to just view
TITLE=""

PY=${PY:-python3}
SCRIPT="./plotters/plot_errors_loglog.py"

# Verify CSVs exist
for f in "${ERRORS_CSV}" "${MAXDL_CSV}" "${MAXDA_CSV}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing CSV: ${f}" >&2
    exit 1
  fi
done

# Ensure output dir if saving
if [[ -n "${OUT}" ]]; then
  mkdir -p "$(dirname "${OUT}")"
fi

# Run
if [[ -z "${OUT}" ]]; then
  "${PY}" "${SCRIPT}" \
    --errors_csv "${ERRORS_CSV}" \
    --maxdL_csv  "${MAXDL_CSV}" \
    --maxdA_csv  "${MAXDA_CSV}" \
    --title "${TITLE}"
else
  "${PY}" "${SCRIPT}" \
    --errors_csv "${ERRORS_CSV}" \
    --maxdL_csv  "${MAXDL_CSV}" \
    --maxdA_csv  "${MAXDA_CSV}" \
    --out "${OUT}" \
    --title "${TITLE}"
fi
