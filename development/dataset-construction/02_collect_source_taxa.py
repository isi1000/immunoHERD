#!/usr/bin/env python3
# Lists the distinct antigen taxa among a host's positives, with counts.
# Feeds step 03 (which proteomes to download) and the K threshold.
# No network, read-only, stdlib only.
#
# taxids can arrive as "a | b" from the mining step; each is split and counted
# separately, so a positive counts once per taxon it belongs to.
#
#   python 02_collect_source_taxa.py --species-dir species/Gallus_gallus
#
# in : <species-dir>/input/metadata.tsv
# out: <species-dir>/negative_pool/source_taxa.tsv
#      (taxid, n_pos, organisms, length_min/median/max)

import argparse
import csv
import os
import statistics
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description="List the antigen taxa among a host's positives")
    ap.add_argument("--species-dir", required=True, help="species/<host> (with input/metadata.tsv)")
    args = ap.parse_args()

    meta = os.path.join(args.species_dir, "input", "metadata.tsv")
    if not os.path.isfile(meta):
        raise SystemExit(f"Not found: {meta}")

    out_dir = os.path.join(args.species_dir, "negative_pool")
    os.makedirs(out_dir, exist_ok=True)
    out_taxa = os.path.join(out_dir, "source_taxa.tsv")

    print(f"Reading: {meta}")

    n_pos_rows = 0
    n_pos_sin_taxid = 0
    # per taxon: number of positives, set of organisms, list of lengths
    tax_npos = defaultdict(int)
    tax_orgs = defaultdict(set)
    tax_lens = defaultdict(list)

    with open(meta, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("class") or "").strip() != "positive":
                continue
            n_pos_rows += 1

            raw_tax = (row.get("taxid") or "").strip()
            taxids = [t.strip() for t in raw_tax.split("|")]
            taxids = [t for t in taxids if t and t != "NA"]

            if not taxids:
                n_pos_sin_taxid += 1
                continue

            organism = (row.get("organism") or "").strip()
            try:
                length = int(float(row.get("length") or ""))
            except (ValueError, TypeError):
                length = None

            for t in taxids:
                tax_npos[t] += 1
                if organism and organism != "NA":
                    tax_orgs[t].add(organism)
                if length is not None:
                    tax_lens[t].append(length)

    # sorted by number of positives, descending
    taxa = sorted(tax_npos.keys(), key=lambda t: tax_npos[t], reverse=True)

    with open(out_taxa, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxid", "n_pos", "organisms",
                    "length_min", "length_median", "length_max"])
        for t in taxa:
            lens = tax_lens[t]
            lmin = min(lens) if lens else "NA"
            lmed = int(statistics.median(lens)) if lens else "NA"
            lmax = max(lens) if lens else "NA"
            orgs = " | ".join(sorted(tax_orgs[t]))
            w.writerow([t, tax_npos[t], orgs, lmin, lmed, lmax])

    # summary
    print("\n--- Antigen taxa summary ---")
    print(f"Positive rows       : {n_pos_rows}")
    print(f"Positives w/o taxid : {n_pos_sin_taxid}")
    print(f"Distinct taxa       : {len(taxa)}")
    for k in (1, 2, 3, 5, 10, 20):
        sel = [t for t in taxa if tax_npos[t] >= k]
        cubren = sum(tax_npos[t] for t in sel)
        print(f"  taxa with >= {k:2d} positives: {len(sel):4d}  (covering {cubren} positives)")

    print("\nTop 15 taxa:")
    print(f"  {'taxid':>10}  {'n_pos':>6}  organism")
    for t in taxa[:15]:
        org = next(iter(sorted(tax_orgs[t])), "")
        print(f"  {t:>10}  {tax_npos[t]:>6}  {org}")

    print(f"\nSaved: {out_taxa}")


if __name__ == "__main__":
    main()
