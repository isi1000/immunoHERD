#!/usr/bin/env bash
#SBATCH --job-name=cluster
#SBATCH --partition=express
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cluster_%j.out
#SBATCH --error=logs/cluster_%j.err
#
# Clusters all curated sequences (positives + negatives) by homology with
# DIAMOND linclust, for the homology-aware CV (no leakage between folds).
#
# approx-id 0 / member-cover 0 groups any pair sharing k-mer seeds. Very
# aggressive on purpose: avoids leakage even below 30% identity.
#
#   bash scripts/06_cluster_sequences.sh --species-dir species/Gallus_gallus
#
# in : <species-dir>/curated/fasta/all.fasta   (from step 05)
# out: <species-dir>/splits/clusters.tsv       (representative <TAB> member)

set -euo pipefail

SPECIES_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --species-dir) SPECIES_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done
[[ -z "$SPECIES_DIR" ]] && { echo "Usage: $0 --species-dir species/<host>"; exit 1; }

command -v diamond >/dev/null 2>&1 || { echo "ERROR: diamond not on PATH"; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-16}"
ALL_FASTA="$SPECIES_DIR/curated/fasta/all.fasta"
SPLITS_DIR="$SPECIES_DIR/splits"
TMP_DIR="$SPECIES_DIR/tmp/clustering"
CLUSTERS="$SPLITS_DIR/clusters.tsv"
SUMMARY="$SPLITS_DIR/clustering_summary.txt"

[[ -f "$ALL_FASTA" ]] || { echo "ERROR: $ALL_FASTA not found"; exit 1; }
mkdir -p "$SPLITS_DIR" "$TMP_DIR"

N=$(grep -c "^>" "$ALL_FASTA" || true)
echo "Sequences to cluster: $N"

diamond linclust \
    --db "$ALL_FASTA" \
    --out "$CLUSTERS" \
    --approx-id 0 \
    --member-cover 0 \
    --memory-limit 30G \
    --threads "$THREADS" \
    --tmpdir "$TMP_DIR"

N_REPS=$(cut -f1 "$CLUSTERS" | sort -u | wc -l)
N_MEMBERS=$(cut -f2 "$CLUSTERS" | sort -u | wc -l)
LARGEST=$(cut -f1 "$CLUSTERS" | sort | uniq -c | sort -rn | head -1 | awk '{print $1}')

cat > "$SUMMARY" <<EOF
Clustering summary
==================
Species dir : $SPECIES_DIR
Input       : $ALL_FASTA ($N sequences)
Clusters    : $CLUSTERS

Representatives (clusters) : $N_REPS
Members                    : $N_MEMBERS
Largest cluster            : $LARGEST
EOF

echo
echo "OK"
echo "Clusters: $CLUSTERS"
cat "$SUMMARY"
