#!/usr/bin/env bash
set -Eeuo pipefail

# Default params (can be overridden via env or CLI)
EXE="./build/euler_feedback"                 # 실행 파일 경로
OUTDIR="./results/trajectory_data"                # 출력 루트 디렉토리
METHODS="all"     # "all" 또는 콤마로 구분: vanilla,feedback,adaptive,SV
HS=("1e-3")                  # 스텝 사이즈 목록: ("1e-3" "1e-4")
BUILD=1                      # 1로 두면 CMake 빌드 수행
VERBOSE=0                    # 1이면 실행 커맨드 echo

usage() {
  cat <<'USAGE'
Usage: ./run_save.sh [--exe ./main] [--outdir traj] [--methods vanilla,SV|all] [--hs 1e-3,1e-4] [--build]
Environment variables also supported: EXE, OUTDIR, METHODS, HS, BUILD=1

Examples:
  ./run_save.sh
  ./run_save.sh --methods all --hs 1e-3,1e-4
  METHODS=feedback,adaptive HS=5e-4 ./run_save.sh
  ./run_save.sh --build                        # build then run with defaults
USAGE
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe) EXE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --methods) METHODS="$2"; shift 2;;
    --hs) HS="$2"; shift 2;;
    --build) BUILD=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

if [[ "${BUILD}" == "1" ]]; then
  if [[ -f CMakeLists.txt ]]; then
    echo "[build] Configuring and building in ./build (Release)..."
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    # Adjust EXE to build output if not explicitly set
    if [[ ! -x "${EXE}" ]]; then
      if [[ -x build/main ]]; then EXE="build/main"; fi
    fi
  else
    echo "[build] No CMakeLists.txt. Skipping build."
  fi
fi

# If EXE not executable, try common locations
if [[ ! -x "${EXE}" ]]; then
  if [[ -x ./main ]]; then EXE="./main"
  elif [[ -x build/main ]]; then EXE="build/main"
  fi
fi

if [[ ! -x "${EXE}" ]]; then
  echo "Error: executable not found at '${EXE}'. Build first or pass --exe."
  exit 2
fi

mkdir -p "${OUTDIR}"

IFS=',' read -ra hlist <<< "${HS}"

for h in "${hlist[@]}"; do
  echo "[run] ${EXE} --save=1 --save_h=${h} --save_methods=${METHODS} --outdir=${OUTDIR}"
  "${EXE}" --save=1 --save_h="${h}" --save_methods="${METHODS}" --outdir="${OUTDIR}"
done

echo "[done] CSVs under '${OUTDIR}/<method>/h=<value>/traj.csv'"
