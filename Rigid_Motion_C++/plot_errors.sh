#!/usr/bin/env bash
set -Eeuo pipefail

# ---- Inputs (produced by your C++ runner) ----
ERRORS_CSV="results/error_data/errors_vs_h.csv"
MAXDE_CSV="results/error_data/maxdE_vs_h.csv"
MAXDPI_CSV="results/error_data/maxdPi_vs_h.csv"
MAXDDET_CSV="results/error_data/maxdDet_vs_h.csv"

# ---- Output (empty => show window) ----
OUT="results/plots/error/error_all.pdf"   # e.g., .png/.pdf; set empty string "" to just view
TITLE=""

PY=${PY:-python3}
SCRIPT="./plotters/plot_errors_loglog.py"

# Ensure output directory if saving
if [[ -n "${OUT}" ]]; then
  mkdir -p "$(dirname "${OUT}")"
fi

# Verify CSVs exist (fail fast if missing)
for f in "${ERRORS_CSV}" "${MAXDE_CSV}" "${MAXDPI_CSV}" "${MAXDDET_CSV}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing CSV: ${f}" >&2
    exit 1
  fi
done

# Run
if [[ -z "${OUT}" ]]; then
  "${PY}" "${SCRIPT}" \
    --errors_csv  "${ERRORS_CSV}" \
    --maxdE_csv   "${MAXDE_CSV}" \
    --maxdPi_csv  "${MAXDPI_CSV}" \
    --maxdDet_csv "${MAXDDET_CSV}" \
    --title "${TITLE}"
else
  "${PY}" "${SCRIPT}" \
    --errors_csv  "${ERRORS_CSV}" \
    --maxdE_csv   "${MAXDE_CSV}" \
    --maxdPi_csv  "${MAXDPI_CSV}" \
    --maxdDet_csv "${MAXDDET_CSV}" \
    --out "${OUT}" \
    --title "${TITLE}"
fi
