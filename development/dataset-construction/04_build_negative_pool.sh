#!/usr/bin/env bash
#SBATCH --job-name=neg_pool
#SBATCH --partition=express
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/neg_pool_%j.out
#SBATCH --error=logs/neg_pool_%j.err
#
# Builds the pool of negative candidates for one host, from the reference
# proteomes downloaded in step 03:
#   1. concatenate the host's proteomes (proteomes_needed.txt)
#   2. drop identical sequences (seqkit rmdup -s)
#   3. drop anything homologous to ANY positive (diamond blastp)
#      -> this also removes candidates that ARE a positive (100% identity)
#
# No network, so it runs either as a SLURM job or straight on a login node.
# Needs diamond and seqkit on the PATH (conda env, or load the modules below).
#
#   bash scripts/04_build_negative_pool.sh --species-dir species/Gallus_gallus
#
# Better: filter homology GLOBALLY by passing a FASTA with the positives of
# ALL hosts. Without --positives it only uses this host's input/pos.fasta.
#   bash scripts/04_build_negative_pool.sh --species-dir species/Gallus_gallus \
#        --positives tmp/all_species_positives.fasta
#
# in : <species-dir>/negative_pool/proteomes_needed.txt   (from step 03)
#      <species-dir>/input/pos.fasta  (or a global --positives)
# out: <species-dir>/negative_pool/candidates.no_homology.fasta

set -euo pipefail

# ── Homology thresholds ─────────────────────────────────────────────────────
EVALUE=0.05
MIN_ID=20          # minimum identity (%)
MIN_QCOV=30        # minimum query coverage (%)

# ── Arguments ───────────────────────────────────────────────────────────────
SPECIES_DIR=""
POSITIVES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --species-dir) SPECIES_DIR="$2"; shift 2 ;;
        --positives)   POSITIVES="$2";   shift 2 ;;
        --evalue)      EVALUE="$2";       shift 2 ;;
        --min-id)      MIN_ID="$2";       shift 2 ;;
        --min-qcov)    MIN_QCOV="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$SPECIES_DIR" ]]; then
    echo "Usage: $0 --species-dir species/<host> [--positives <global_pos.fasta>]"
    exit 1
fi

NEEDED="$SPECIES_DIR/negative_pool/proteomes_needed.txt"
OUTDIR="$SPECIES_DIR/negative_pool"
[[ -z "$POSITIVES" ]] && POSITIVES="$SPECIES_DIR/input/pos.fasta"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ── Modules (uncomment if not using a conda env) ────────────────────────────
# module load DIAMOND
# module load SeqKit

for tool in diamond seqkit; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: '$tool' not on PATH. Activate the environment or load the module."
        exit 1
    fi
done

for f in "$NEEDED" "$POSITIVES"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f"; exit 1; }
done

TMPDIR="$SPECIES_DIR/tmp/negative_pool"
mkdir -p "$OUTDIR" "$TMPDIR"

RAW="$TMPDIR/pool.raw.fasta"
CLEAN="$TMPDIR/pool.clean.fasta"
DB="$TMPDIR/positives_db"
HITS="$TMPDIR/hits.tsv"
HOMOLOGS="$TMPDIR/homologous_ids.txt"
OUT_FASTA="$OUTDIR/candidates.no_homology.fasta"
SUMMARY="$OUTDIR/pool_summary.txt"

echo "========================================"
echo "  Building the negative pool"
echo "  needed     : $NEEDED"
echo "  positives  : $POSITIVES"
echo "  out-dir    : $OUTDIR"
echo "  threads    : $THREADS"
echo "  homology   : evalue<=$EVALUE id>=$MIN_ID qcov>=$MIN_QCOV"
echo "========================================"

# ── 1. Concatenate the host proteomes ───────────────────────────────────────
echo "[1/5] Concatenating proteomes..."
: > "$RAW"
n_files=0
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    if [[ -s "$f" ]]; then
        zcat "$f" >> "$RAW"
        n_files=$((n_files + 1))
    else
        echo "  warning: missing or empty: $f" >&2
    fi
done < "$NEEDED"
N_RAW=$(grep -c "^>" "$RAW" || true)
echo "  proteomes joined: $n_files | sequences: $N_RAW"

# ── 2. Drop exact duplicates ────────────────────────────────────────────────
echo "[2/5] Dropping exact duplicates (rmdup -s)..."
seqkit rmdup -s "$RAW" -o "$CLEAN" 2> "$TMPDIR/rmdup.log" || true
N_CLEAN=$(grep -c "^>" "$CLEAN" || true)
echo "  after dedup: $N_CLEAN"

# ── 3. DIAMOND database of positives ────────────────────────────────────────
echo "[3/5] Building the DIAMOND database of positives..."
diamond makedb --in "$POSITIVES" -d "$DB" >/dev/null 2> "$TMPDIR/makedb.log"
N_POS=$(grep -c "^>" "$POSITIVES" || true)
echo "  positives in the database: $N_POS"

# ── 4. Search for homologs of the positives ─────────────────────────────────
echo "[4/5] DIAMOND blastp (pool vs positives)..."
diamond blastp \
    -q "$CLEAN" \
    -d "$DB" \
    -o "$HITS" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp \
    --evalue "$EVALUE" \
    --id "$MIN_ID" \
    --query-cover "$MIN_QCOV" \
    --max-target-seqs 1 \
    --threads "$THREADS" \
    >/dev/null 2> "$TMPDIR/blastp.log"

cut -f1 "$HITS" | sort -u > "$HOMOLOGS"
N_HOMOLOGS=$(wc -l < "$HOMOLOGS")
echo "  candidates homologous to a positive: $N_HOMOLOGS"

# ── 5. Remove homologs -> final pool ────────────────────────────────────────
echo "[5/5] Removing homologs..."
if [[ -s "$HOMOLOGS" ]]; then
    seqkit grep -v -f "$HOMOLOGS" "$CLEAN" -o "$OUT_FASTA" 2> "$TMPDIR/grep.log"
else
    cp "$CLEAN" "$OUT_FASTA"
fi
N_FINAL=$(grep -c "^>" "$OUT_FASTA" || true)

# ── Summary ─────────────────────────────────────────────────────────────────
cat > "$SUMMARY" <<EOF
Negative pool summary
=====================
Needed list : $NEEDED
Positives   : $POSITIVES ($N_POS sequences)
Homology    : evalue<=$EVALUE id>=$MIN_ID qcov>=$MIN_QCOV

Proteomes joined         : $n_files
Raw sequences            : $N_RAW
After exact dedup        : $N_CLEAN
Homologous to positives  : $N_HOMOLOGS
Final candidate pool     : $N_FINAL

Output: $OUT_FASTA
EOF

echo
echo "OK"
echo "Pool: $OUT_FASTA ($N_FINAL sequences)"
echo "Summary: $SUMMARY"
