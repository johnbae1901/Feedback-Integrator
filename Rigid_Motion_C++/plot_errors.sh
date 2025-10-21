#!/usr/bin/env bash
set -Eeuo pipefail

CSV="error_data/errors_vs_h.csv"
OUT=""            # leave empty to view only; set like OUT="plots/error_curve.png" or ".pdf"
TITLE="Error vs Step Size (log-log)"

PY="python3"
SCRIPT="./plotters/plot_errors_loglog.py"

mkdir -p "$(dirname "${CSV}")" || true

if [[ -z "${OUT}" ]]; then
  "${PY}" "${SCRIPT}" --csv "${CSV}" --title "${TITLE}"
else
  mkdir -p "$(dirname "${OUT}")"
  "${PY}" "${SCRIPT}" --csv "${CSV}" --out "${OUT}" --title "${TITLE}"
fi
