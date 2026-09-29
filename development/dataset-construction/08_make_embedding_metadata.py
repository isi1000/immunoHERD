#!/usr/bin/env python3
# Builds the input for the ESM2 embeddings (run in Colab): one FASTA with ALL
# curated sequences and a master metadata table with both CV columns attached
# (homology-aware and leave-one-pathogen-out).
#
# No network, read-only, stdlib only.
#
#   python 08_make_embedding_metadata.py --species-dir species/Gallus_gallus
#
# in : <sd>/curated/sequences.tsv
#      <sd>/splits/cv_homology.tsv   (from step 07)
#      <sd>/splits/cv_lopo.tsv       (from step 07)
# out: <sd>/embeddings/esm2/all_sequences.fasta
#      <sd>/embeddings/esm2/all_embedding_metadata.tsv

import argparse
import csv
import os
from collections import Counter


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description="Metadata + FASTA for the ESM2 embeddings")
    ap.add_argument("--species-dir", required=True)
    args = ap.parse_args()

    sd = args.species_dir
    seq_file = os.path.join(sd, "curated", "sequences.tsv")
    homo_file = os.path.join(sd, "splits", "cv_homology.tsv")
    lopo_file = os.path.join(sd, "splits", "cv_lopo.tsv")
    out_dir = os.path.join(sd, "embeddings", "esm2")
    os.makedirs(out_dir, exist_ok=True)
    out_meta = os.path.join(out_dir, "all_embedding_metadata.tsv")
    out_fasta = os.path.join(out_dir, "all_sequences.fasta")

    for f in (seq_file, homo_file, lopo_file):
        if not os.path.isfile(f):
            raise SystemExit(f"Not found: {f} (did you run 05/06/07?)")

    records = read_tsv(seq_file)

    homo = {r["id"]: r for r in read_tsv(homo_file)}
    lopo = {r["id"]: r for r in read_tsv(lopo_file)}

    cols = ["id", "class", "label", "confidence_level", "type_group",
            "organism", "taxid", "length", "length_bin",
            "cv_homology_fold", "cluster_id",
            "cv_lopo_group", "genus", "source"]

    out_rows = []
    for r in records:
        sid = r["id"]
        h = homo.get(sid, {})
        l = lopo.get(sid, {})
        out_rows.append({
            "id": sid,
            "class": r.get("class", "NA"),
            "label": r.get("label", "NA"),
            "confidence_level": r.get("confidence_level", "NA"),
            "type_group": r.get("type_group", "NA"),
            "organism": r.get("organism", "NA"),
            "taxid": r.get("taxid", "NA"),
            "length": r.get("length", "NA"),
            "length_bin": r.get("length_bin", "NA"),
            "cv_homology_fold": h.get("fold", "NA"),
            "cluster_id": h.get("cluster_id", "NA"),
            "cv_lopo_group": l.get("group", "NA"),
            "genus": l.get("genus", "NA"),
            "source": r.get("source", "NA"),
        })

    with open(out_meta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # FASTA in the same order as the metadata
    n_fasta = 0
    with open(out_fasta, "w", encoding="utf-8") as f:
        for r in records:
            seq = r.get("sequence", "")
            if not seq or seq == "NA":
                continue
            f.write(f">{r['id']}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
            n_fasta += 1

    print(f"[OK] sequences: {len(records)} | FASTA written: {n_fasta}")
    print(f"  {out_meta}")
    print(f"  {out_fasta}")
    print("\nBy class:", dict(Counter(r["class"] for r in out_rows)))
    print("Homology folds:", dict(Counter(r["cv_homology_fold"] for r in out_rows)))
    print("LOPO groups:", len(set(r["cv_lopo_group"] for r in out_rows)))


if __name__ == "__main__":
    main()
