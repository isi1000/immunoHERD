#!/usr/bin/env bash
# Wrapper over run_pipeline.sh: runs it for EVERY host under species/*,
# reusing a single global FASTA of positives instead of letting each host
# rebuild its own, and without aborting the whole batch when one host fails
# (typically at step 07).
#
#   bash scripts/run_all_species.sh
#   bash scripts/run_all_species.sh -r 10
#   bash scripts/run_all_species.sh -r 10 --skip-download
#
# Hosts run one after another, not in parallel, because the resources
# available on the cluster for this job are not known in advance. The
# expensive shared part (the global positives) is already computed once, so
# switching to background jobs or a SLURM array only needs the loop changed.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# ── Defaults (same names and flags as run_pipeline.sh) ──────────────────────
RATIO=1
MIN_POS=1
FOLDS=5
MIN_GROUP=10
SKIP_DOWNLOAD=0
RUN_BIAS=1

usage() {
    cat <<EOF
Usage: bash scripts/run_all_species.sh [options]
  -r, --ratio       Negatives per positive (step 05)                 [def: $RATIO]
  -k, --min-pos     K: download taxa with >=K positives (03)         [def: $MIN_POS]
  -f, --folds       Folds of the homology-aware CV (07)              [def: $FOLDS]
  -g, --min-group   Minimum positives per LOPO group (07)            [def: $MIN_GROUP]
      --skip-download   Skip step 03 (use the existing proteome cache)
      --no-bias         Skip the bias QC (step 09)
  -h, --help        This help

The same parameters are applied to every host found under species/*/.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ratio)      RATIO="$2"; shift 2 ;;
        -k|--min-pos)    MIN_POS="$2"; shift 2 ;;
        -f|--folds)      FOLDS="$2"; shift 2 ;;
        -g|--min-group)  MIN_GROUP="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        --no-bias)       RUN_BIAS=0; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

# ── 1. Global positives FASTA: built ONCE for every host ────────────────────
#    run_pipeline.sh rebuilds it on every call; here it is fixed once and
#    passed to each host with --positives, so the work is not repeated N times
#    and concurrent runs cannot write over the same tmp file.
GLOBAL_POS="tmp/all_species_positives.fasta"
mkdir -p tmp
echo "============================================================"
echo ">>> Building the global positives FASTA (once) -> $GLOBAL_POS"
: > "$GLOBAL_POS"
n_sp=0
for pf in species/*/input/pos.fasta; do
    [[ -f "$pf" ]] || continue
    cat "$pf" >> "$GLOBAL_POS"
    n_sp=$((n_sp + 1))
    echo "    + $pf"
done
if [[ "$n_sp" -eq 0 ]]; then
    echo "ERROR: no species/*/input/pos.fasta found"; exit 1
fi
seqkit rmdup -n "$GLOBAL_POS" -o "$GLOBAL_POS.dedup" 2>/dev/null && mv "$GLOBAL_POS.dedup" "$GLOBAL_POS"
N_GP=$(grep -c "^>" "$GLOBAL_POS" || true)
echo "    hosts with pos.fasta: $n_sp | unique positives (by ID): $N_GP"
echo "============================================================"

# ── 2. Loop over every host whose input is ready (steps 00/01) ──────────────
FAILED=()
OK=()
SKIPPED=()

for dir in species/*/; do
    sp="$(basename "$dir")"
    META="$dir/input/metadata.tsv"
    POS="$dir/input/pos.fasta"
    if [[ ! -f "$META" || ! -f "$POS" ]]; then
        echo ">>> SKIP $sp (no input/metadata.tsv or input/pos.fasta - run 00/01 first)"
        SKIPPED+=("$sp")
        continue
    fi

    echo
    echo "################################################################"
    echo "### HOST: $sp"
    echo "################################################################"

    ARGS=(-s "$sp" -r "$RATIO" -k "$MIN_POS" -f "$FOLDS" -g "$MIN_GROUP" -p "$GLOBAL_POS")
    [[ "$SKIP_DOWNLOAD" -eq 1 ]] && ARGS+=(--skip-download)
    [[ "$RUN_BIAS" -eq 0 ]] && ARGS+=(--no-bias)

    # run_pipeline.sh has its own 'set -e': a mid-run failure (e.g. step 07
    # with a LOPO group below --min-group) must not kill the whole loop, so it
    # is caught here and the batch carries on with the next host.
    if bash scripts/run_pipeline.sh "${ARGS[@]}"; then
        OK+=("$sp")
    else
        echo ">>> FAILED on $sp (check the log above, especially step 07) - carrying on"
        FAILED+=("$sp")
    fi
done

# ── 3. Final summary ────────────────────────────────────────────────────────
echo
echo "============================================================"
echo "  run_all_species.sh summary"
echo "============================================================"
echo "Completed OK (${#OK[@]}): ${OK[*]:-none}"
echo "Failed       (${#FAILED[@]}): ${FAILED[*]:-none}"
echo "Skipped      (${#SKIPPED[@]}): ${SKIPPED[*]:-none}"
echo "============================================================"

[[ "${#FAILED[@]}" -eq 0 ]]
