#!/usr/bin/env bash
# Runs the negatives + splits + embeddings pipeline for ONE host, chaining
# steps 02 -> 09 with a chosen neg:pos ratio.
#
# Requires the positives to exist already (steps 00 and 01):
#   species/<HOST>/input/metadata.tsv
#   species/<HOST>/input/pos.fasta
#
# Needs python (+ete3), diamond and seqkit on the PATH. Only step 03 needs
# INTERNET; use --skip-download to run 04->09 on a node without network, once
# the proteome cache is in place.
#
#   bash scripts/run_pipeline.sh -s Gallus_gallus -r 10
#   bash scripts/run_pipeline.sh --species Gallus_gallus --ratio 1 --min-pos 1
#   bash scripts/run_pipeline.sh -s Gallus_gallus -r 10 --positives tmp/all_species_positives.fasta
#   bash scripts/run_pipeline.sh -s Gallus_gallus -r 10 --skip-download

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
SPECIES=""
RATIO=1            # negatives per positive (step 05)
MIN_POS=1          # K: download taxa with >=K positives (step 03)
FOLDS=5            # folds of the homology-aware CV (step 07)
MIN_GROUP=10       # minimum positives per LOPO group (step 07)
POSITIVES=""       # global positives FASTA for the homology filter (step 04);
                   # empty = built by concatenating every host's positives
SKIP_DOWNLOAD=0    # 1 = skip step 03 (use the cache already downloaded)
RUN_BIAS=1         # 1 = run 09_check_bias

# ── Argument parsing ────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: bash scripts/run_pipeline.sh -s <HOST> [options]

  -s, --species     Host folder under species/ (e.g. Gallus_gallus)      [required]
  -r, --ratio       Negatives per positive (step 05)                     [def: $RATIO]
  -k, --min-pos     K: download taxa with >=K positives (step 03)        [def: $MIN_POS]
  -f, --folds       Folds of the homology-aware CV (step 07)             [def: $FOLDS]
  -g, --min-group   Minimum positives per LOPO group (step 07)           [def: $MIN_GROUP]
  -p, --positives   Global positives FASTA for homology (step 04)        [def: auto = every host]
      --skip-download   Skip step 03 (use the existing proteome cache)
      --no-bias         Skip the bias QC (step 09)
  -h, --help        This help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--species)    SPECIES="$2"; shift 2 ;;
        -r|--ratio)      RATIO="$2"; shift 2 ;;
        -k|--min-pos)    MIN_POS="$2"; shift 2 ;;
        -f|--folds)      FOLDS="$2"; shift 2 ;;
        -g|--min-group)  MIN_GROUP="$2"; shift 2 ;;
        -p|--positives)  POSITIVES="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        --no-bias)       RUN_BIAS=0; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$SPECIES" ]]; then
    echo "ERROR: --species is required"; usage; exit 1
fi

# ── Move to the project root (parent of scripts/) ───────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

SPECIES_DIR="species/$SPECIES"
META="$SPECIES_DIR/input/metadata.tsv"
POS="$SPECIES_DIR/input/pos.fasta"

# ── Pre-flight checks ───────────────────────────────────────────────────────
echo "============================================================"
echo "  PIPELINE — $SPECIES"
echo "  project root  : $PROJECT_ROOT"
echo "  ratio neg:pos : 1:$RATIO   | min-pos(K): $MIN_POS"
echo "  folds: $FOLDS | min-group LOPO: $MIN_GROUP"
echo "  homology positives: ${POSITIVES:-auto (every host)}"
echo "  skip-download: $SKIP_DOWNLOAD | bias QC: $RUN_BIAS"
echo "============================================================"

for f in "$META" "$POS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f (did you run 00/01?)"; exit 1; }
done
for t in python diamond seqkit; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not on PATH"; exit 1; }
done

# ── Helper: run one step with a banner and a timer ──────────────────────────
run_step() {
    local title="$1"; shift
    echo
    echo "------------------------------------------------------------"
    echo ">>> $title"
    echo "    \$ $*"
    echo "------------------------------------------------------------"
    local t0=$SECONDS
    "$@"
    echo "    [OK] $title  (${SECONDS}s total, +$((SECONDS - t0))s)"
}

# ── Pipeline ──────────────────────────────────────────────────────────────────
run_step "02 - Collect antigen taxa" \
    python scripts/02_collect_source_taxa.py --species-dir "$SPECIES_DIR"

if [[ "$SKIP_DOWNLOAD" -eq 0 ]]; then
    run_step "03 - Resolve and download reference proteomes (NEEDS INTERNET)" \
        python scripts/03_download_proteomes.py --species-dir "$SPECIES_DIR" --min-pos "$MIN_POS"
else
    echo; echo ">>> 03 - SKIPPED (--skip-download): using the existing databases/proteomes cache"
fi

# ── Shared positives FASTA across ALL hosts (cross-host homology) ───────────
# Without --positives it is built by concatenating species/*/input/pos.fasta,
# so step 04 drops from the negative pool anything homologous to a positive of
# ANY host, not just this one.
if [[ -z "$POSITIVES" ]]; then
    GLOBAL_POS="tmp/all_species_positives.fasta"
    mkdir -p tmp
    echo
    echo ">>> Building the shared positives FASTA (all hosts) -> $GLOBAL_POS"
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
    # Drop duplicate IDs (the same epitope present in several hosts), or
    # 'diamond makedb' fails on repeated accessions. Keeps the first one.
    seqkit rmdup -n "$GLOBAL_POS" -o "$GLOBAL_POS.dedup" 2>/dev/null && mv "$GLOBAL_POS.dedup" "$GLOBAL_POS"
    N_GP=$(grep -c "^>" "$GLOBAL_POS" || true)
    echo "    hosts with pos.fasta: $n_sp | unique positives (by ID): $N_GP"
    POSITIVES="$GLOBAL_POS"
fi

POS_ARGS=()
[[ -n "$POSITIVES" ]] && POS_ARGS=(--positives "$POSITIVES")
run_step "04 - Build the negative pool (dedup + homology)" \
    bash scripts/04_build_negative_pool.sh --species-dir "$SPECIES_DIR" "${POS_ARGS[@]}"

run_step "05 - Match negatives (ratio 1:$RATIO)" \
    python scripts/05_match_negatives.py --species-dir "$SPECIES_DIR" --ratio "$RATIO"

run_step "06 - Homology clustering (DIAMOND linclust)" \
    bash scripts/06_cluster_sequences.sh --species-dir "$SPECIES_DIR"

run_step "07 - CV splits (homology-aware + LOPO)" \
    python scripts/07_make_cv_splits.py --species-dir "$SPECIES_DIR" --folds "$FOLDS" --min-group "$MIN_GROUP"

run_step "08 - Metadata + FASTA for the ESM2 embeddings" \
    python scripts/08_make_embedding_metadata.py --species-dir "$SPECIES_DIR"

if [[ "$RUN_BIAS" -eq 1 ]]; then
    run_step "09 - Bias audit (QC)" \
        python scripts/09_check_bias.py --species-dir "$SPECIES_DIR"
fi

# ── Summary ─────────────────────────────────────────────────────────────────
echo
echo "============================================================"
echo "  PIPELINE DONE — $SPECIES (ratio 1:$RATIO) in ${SECONDS}s"
echo "============================================================"
echo "Main outputs:"
echo "  $SPECIES_DIR/negative_pool/candidates.no_homology.fasta"
echo "  $SPECIES_DIR/negative_pool/matched_negatives.fasta + matching.tsv"
echo "  $SPECIES_DIR/curated/sequences.tsv + curated/fasta/all.fasta"
echo "  $SPECIES_DIR/splits/cv_homology.tsv + cv_lopo.tsv"
echo "  $SPECIES_DIR/embeddings/esm2/all_sequences.fasta + all_embedding_metadata.tsv"
[[ "$RUN_BIAS" -eq 1 ]] && echo "  $SPECIES_DIR/qc/bias_report.txt"
