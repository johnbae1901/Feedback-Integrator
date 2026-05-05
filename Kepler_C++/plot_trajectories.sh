#!/usr/bin/env bash
set -Eeuo pipefail
PY="python3"
PYSCRIPT="./plotters/plot_trajectories.py"
BASE="results/trajectory_data"
METHODS=("vanilla" "vanilla_feedback" "feedback" "adaptive" "SV")
# METHODS=("SV")
HS=("1e-3")
OUTROOT="results/plots/trajectories"
OUTNAME="traj_x1x2"
RASTERIZE_TRAJ="${RASTERIZE_TRAJ:-10}"
RASTER_DPI="${RASTER_DPI:-600}"
for method in "${METHODS[@]}"; do
  for h in "${HS[@]}"; do
    outdir="${OUTROOT}/${method}/h=${h}"
    mkdir -p "${outdir}"
    outpdf="${outdir}/${OUTNAME}.pdf"
    "${PY}" "${PYSCRIPT}" --base "${BASE}" --h "${h}" --method "${method}" --out "${outpdf}" --rasterize_traj "${RASTERIZE_TRAJ}" --raster_dpi "${RASTER_DPI}"
  done
done
