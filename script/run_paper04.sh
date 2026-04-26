#!/usr/bin/env bash
# run_paper04.sh
# Run all paper04 scripts (proofs/ and paper04/).
# Each script writes its own timestamped log to the paper-level data/ directory.
#
# Usage (from the script/ directory):
#   bash run_paper04.sh                          # run all -> ../data/
#   bash run_paper04.sh proofs                   # run proofs/ only
#   bash run_paper04.sh paper04                  # run paper04/ only
#   bash run_paper04.sh --output-dir my_logs     # custom output dir
#   bash run_paper04.sh -o my_logs paper04       # combined
#
# Output: <output-dir>/<script_name>_<timestamp>.log  (written by each script)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# ── argument parsing ──────────────────────────────────────────────────────────

TARGET="all"
DATA_DIR="${SCRIPT_DIR}/../data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir|-o)
            DATA_DIR="$2"
            shift 2
            ;;
        all|proofs|paper04)
            TARGET="$1"
            shift
            ;;
        *)
            echo "Usage: bash run_paper04.sh [--output-dir DIR] [all|proofs|paper04]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

export PYTHONIOENCODING=utf-8
export PYTHONUTF8=1
export DPPU_LOG_DIR="${DATA_DIR}"
export DPPU_LOG_STDOUT=0

# ── helper ────────────────────────────────────────────────────────────────────

FAILED_SCRIPTS=()

run_script() {
    local rel_path="$1"
    local script_file="${SCRIPTS_DIR}/${rel_path}"
    local base
    base="$(basename "${rel_path}" .py)"
    echo ">>> ${rel_path}"
    local t0
    t0=$(date +%s)

    if python "${script_file}"; then
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    OK  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
    else
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    FAILED  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
        FAILED_SCRIPTS+=("${rel_path}")
    fi
}

# ── script lists ──────────────────────────────────────────────────────────────

# Proof scripts: symbolic / algebraic proofs (fast)
PROOFS_SCRIPTS=(
    proofs/sol3_structure.py
    proofs/cs_cancellation.py
    proofs/kk_higgsing.py
    proofs/eta_kinetic_from_contortion.py
    proofs/weyl_scalar.py
    proofs/aps_zero_t3_s3.py
    proofs/landau_levels_nil3.py
    proofs/eta_aps_sol3.py
    proofs/kk_normalization.py
)

# paper04 scripts: defect localization, APS, EC minima (heavier computation)
PAPER04_SCRIPTS=(
    paper04/ec_slice_minima.py
    paper04/eta_aps_nil3.py
    paper04/eta_defect_coefficients.py
    paper04/torsional_charge.py
    paper04/defect_localization.py
)

# ── dispatch ──────────────────────────────────────────────────────────────────

echo "========================================"
echo "  run_paper04.sh  (target: ${TARGET})"
echo "  output dir: ${DATA_DIR}"
echo "========================================"
echo ""

if [[ "${TARGET}" == "all" || "${TARGET}" == "proofs" ]]; then
    echo "--- proofs/ ---"
    for s in "${PROOFS_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

if [[ "${TARGET}" == "all" || "${TARGET}" == "paper04" ]]; then
    echo "--- paper04/ ---"
    for s in "${PAPER04_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

# ── summary ───────────────────────────────────────────────────────────────────

echo "========================================"
if [[ ${#FAILED_SCRIPTS[@]} -eq 0 ]]; then
    echo "  All scripts completed successfully."
else
    echo "  ${#FAILED_SCRIPTS[@]} script(s) failed:"
    for s in "${FAILED_SCRIPTS[@]}"; do
        echo "    - ${s}"
    done
fi
echo "========================================"
