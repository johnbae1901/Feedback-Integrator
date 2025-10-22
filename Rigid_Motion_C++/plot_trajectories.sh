#!/usr/bin/env bash
set -Eeuo pipefail

########################################
# VIEW-ONLY DEFAULTS (no PNG saved)
PY="python3"
PYSCRIPT="./plotters/plot_trajectories.py"
BASE="results/trajectory_data"
METHODS=("vanilla" "vanilla_feedback" "feedback" "adaptive" "strang")
HS=("1e-4")
VERBOSE=0

# If you later want to save, flip this to 1 and set OUTROOT/OUTNAME.
SAVE=1
OUTROOT="results/plots/trajectories"
OUTNAME="traj"
USE_TITLE=0
TITLE_PREFIX="Projection"
########################################

[[ "${VERBOSE}" == "1" ]] && set -x

for method in "${METHODS[@]}"; do
  for h in "${HS[@]}"; do
    if [[ "${SAVE}" == "1" ]]; then
      outdir="${OUTROOT}/${method}/h=${h}"
      mkdir -p "${outdir}"
      outpdf="${outdir}/${OUTNAME}.pdf"
      echo "[plot/save] ${method} h=${h} -> ${outpdf}"
      if [[ "${USE_TITLE}" == "1" ]]; then
        "${PY}" "${PYSCRIPT}" --base "${BASE}" --h "${h}" --method "${method}" --out "${outpdf}" --title "${TITLE_PREFIX} ${method} (h=${h})"
      else
        "${PY}" "${PYSCRIPT}" --base "${BASE}" --h "${h}" --method "${method}" --out "${outpdf}"
      fi
    else
      echo "[plot(view)] ${method} h=${h}"
      if [[ "${USE_TITLE}" == "1" ]]; then
        "${PY}" "${PYSCRIPT}" --base "${BASE}" --h "${h}" --method "${method}" --title "${TITLE_PREFIX} ${method} (h=${h})"
      else
        "${PY}" "${PYSCRIPT}" --base "${BASE}" --h "${h}" --method "${method}"
      fi
    fi
  done
done
