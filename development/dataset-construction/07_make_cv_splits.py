#!/usr/bin/env python3
# Builds TWO cross-validation schemes over the curated dataset (positives +
# matched negatives). There is no held-out "evidence" test set: every fold
# takes its turn as the test fold.
#
#   1) cv_homology.tsv  k-fold over the homology CLUSTERS from step 06,
#      stratified by class and length. No homology leakage between folds.
#
#   2) cv_lopo.tsv      leave-one-pathogen-out at GENUS level. Each negative
#      inherits the genus of the positive it was matched to (matching.tsv).
#      Genera with fewer than MIN_GROUP positives roll up to family -> order
#      -> OTHER_RARE. Each 'group' is a fold and is left out whole.
#
# Needs ete3 (NCBI taxonomy).
#
#   python 07_make_cv_splits.py --species-dir species/Gallus_gallus
#
# in : <sd>/curated/sequences.tsv
#      <sd>/splits/clusters.tsv             (step 06)
#      <sd>/negative_pool/matching.tsv      (step 05)
#      <sd>/negative_pool/proteome_plan.tsv (step 03, taxid remap)
# out: <sd>/splits/cv_homology.tsv
#      <sd>/splits/cv_lopo.tsv
#      <sd>/splits/cv_summary.txt

import argparse
import csv
import os
import random
from collections import Counter, defaultdict

from ete3 import NCBITaxa


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def first_taxid(raw):
    if not raw:
        return None
    for t in str(raw).replace("|", " ").split():
        if t and t != "NA":
            return to_int(t)
    return None


def length_bin(n):
    n = to_int(n) or 0
    if n <= 15:  return "L01_very_short"
    if n <= 30:  return "L02_short_peptide"
    if n <= 60:  return "L03_long_peptide"
    if n <= 150: return "L04_short_protein"
    if n <= 400: return "L05_medium_protein"
    return "L06_long_protein"


def main():
    ap = argparse.ArgumentParser(description="Build the homology-aware and leave-one-pathogen-out CV splits")
    ap.add_argument("--species-dir", required=True)
    ap.add_argument("--folds", type=int, default=5, help="K for the homology-aware CV")
    ap.add_argument("--min-group", type=int, default=10, help="minimum POSITIVES per LOPO group")
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()
    random.seed(args.seed)

    sd = args.species_dir
    seq_file = os.path.join(sd, "curated", "sequences.tsv")
    clusters_file = os.path.join(sd, "splits", "clusters.tsv")
    matching_file = os.path.join(sd, "negative_pool", "matching.tsv")
    plan_file = os.path.join(sd, "negative_pool", "proteome_plan.tsv")
    splits_dir = os.path.join(sd, "splits")
    os.makedirs(splits_dir, exist_ok=True)

    for f in (seq_file, clusters_file):
        if not os.path.isfile(f):
            raise SystemExit(f"Not found: {f}")

    # taxid -> real taxid, from the proteome plan
    remap = {}
    if os.path.isfile(plan_file):
        for row in read_tsv(plan_file):
            if row.get("status") in ("resolved", "uniprotkb_taxon"):
                pt, ot = to_int(row.get("proteome_taxid")), to_int(row.get("taxid"))
                if pt and ot:
                    remap[ot] = pt

    # neg_id -> pos_id, from the matching
    neg2pos = {}
    if os.path.isfile(matching_file):
        for row in read_tsv(matching_file):
            nid = row.get("negative_id")
            if nid and nid != "NA":
                neg2pos[nid] = row.get("positive_id")

    # clusters: member -> representative
    rep_of = {}
    with open(clusters_file, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                rep_of[parts[1]] = parts[0]

    # sequences
    records = read_tsv(seq_file)
    print(f"[1/4] sequences: {len(records)}")

    # ── lineage (genus/family/order) by real taxid ──────────────────────────
    ncbi = NCBITaxa()
    chain_cache = {}

    def chain_for(taxid):
        if taxid is None:
            return {"genus": None, "family": None, "order": None}
        if taxid in chain_cache:
            return chain_cache[taxid]
        out = {"genus": None, "family": None, "order": None}
        try:
            lin = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lin)
            names = ncbi.get_taxid_translator(lin)
            for t in lin:
                r = ranks.get(t)
                if r in out:
                    out[r] = names.get(t, str(t))
        except Exception:  # noqa: BLE001
            pass
        chain_cache[taxid] = out
        return out

    # lineage per POSITIVE, via its real taxid
    pos_chain = {}
    for r in records:
        if r.get("class") == "positive":
            ot = first_taxid(r.get("taxid"))
            real = remap.get(ot, ot)
            pos_chain[r["id"]] = chain_for(real)

    # lineage per sequence; negatives inherit their positive's
    chains = {}
    for r in records:
        sid = r["id"]
        if r.get("class") == "positive":
            chains[sid] = pos_chain.get(sid, {"genus": None, "family": None, "order": None})
        else:
            pid = neg2pos.get(sid)
            if pid and pid in pos_chain:
                chains[sid] = pos_chain[pid]
            else:
                chains[sid] = chain_for(remap.get(first_taxid(r.get("taxid")),
                                                  first_taxid(r.get("taxid"))))

    # ── 2) LOPO by genus, rolling up rare groups (threshold on POSITIVES) ───
    print("[2/4] LOPO by genus (rolling up rare groups by positive count)")
    ids = [r["id"] for r in records]
    is_pos = {r["id"]: (r.get("class") == "positive") for r in records}

    def pos_sizes(grp):
        c = Counter()
        for sid in ids:
            if is_pos[sid] and grp[sid] is not None:
                c[grp[sid]] += 1
        return c

    group = {}
    for sid in ids:
        g = chains[sid]["genus"]
        group[sid] = f"genus:{g}" if g else None

    # promote groups with < min_group POSITIVES: genus -> family -> order
    for level in ("family", "order"):
        sizes = pos_sizes(group)
        labels = {v for v in group.values() if v is not None}
        small = {lab for lab in labels if sizes.get(lab, 0) < args.min_group}
        for sid in ids:
            if group[sid] is None or group[sid] in small:
                v = chains[sid][level]
                if v:
                    group[sid] = f"{level}:{v}"
    sizes = pos_sizes(group)
    labels = {v for v in group.values() if v is not None}
    small = {lab for lab in labels if sizes.get(lab, 0) < args.min_group}
    for sid in ids:
        if group[sid] is None or group[sid] in small:
            group[sid] = "OTHER_RARE"

    n_lopo_groups = len(set(group.values()))
    print(f"  LOPO groups (folds): {n_lopo_groups}")

    by_id = {r["id"]: r for r in records}
    with open(os.path.join(splits_dir, "cv_lopo.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "class", "taxid", "genus", "family", "order", "group"])
        for sid in sorted(ids):
            r = by_id[sid]
            ch = chains[sid]
            w.writerow([sid, r.get("class", "NA"), r.get("taxid", "NA"),
                        ch["genus"] or "NA", ch["family"] or "NA", ch["order"] or "NA",
                        group[sid]])

    # ── 1) homology-aware k-fold over clusters ──────────────────────────────
    print(f"[3/4] homology-aware CV ({args.folds} folds over clusters)")
    for r in records:
        r["_lbin"] = length_bin(r.get("length"))
        r["_rep"] = rep_of.get(r["id"], r["id"])

    clusters = defaultdict(list)
    for r in records:
        clusters[r["_rep"]].append(r)

    cstats = []
    for rep, members in clusters.items():
        cl = Counter(m["class"] for m in members)
        lb = Counter(m["_lbin"] for m in members)
        cstats.append({"rep": rep, "members": members, "size": len(members),
                       "n_pos": cl.get("positive", 0), "n_neg": cl.get("negative", 0),
                       "lbin": lb.most_common(1)[0][0]})
    random.shuffle(cstats)
    cstats.sort(key=lambda x: x["size"], reverse=True)

    total = len(records)
    tot_pos = sum(1 for r in records if r["class"] == "positive")
    tot_neg = total - tot_pos
    t_total, t_pos, t_neg = total / args.folds, max(tot_pos, 1) / args.folds, max(tot_neg, 1) / args.folds
    tot_lb = Counter(r["_lbin"] for r in records)
    t_lb = {k: max(v / args.folds, 1) for k, v in tot_lb.items()}

    folds = [{"name": f"fold_{i+1}", "n": 0, "pos": 0, "neg": 0, "lb": Counter()}
             for i in range(args.folds)]

    def score(fold, c):
        s = fold["n"] / max(t_total, 1)
        if c["n_pos"]:
            s += 3 * fold["pos"] / t_pos
        if c["n_neg"]:
            s += 3 * fold["neg"] / t_neg
        s += fold["lb"][c["lbin"]] / t_lb.get(c["lbin"], 1)
        return s

    assign = {}
    for c in cstats:
        fb = min(folds, key=lambda fo: score(fo, c))
        fb["n"] += c["size"]; fb["pos"] += c["n_pos"]; fb["neg"] += c["n_neg"]
        fb["lb"][c["lbin"]] += c["size"]
        assign[c["rep"]] = fb["name"]

    with open(os.path.join(splits_dir, "cv_homology.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "class", "length_bin", "cluster_id", "fold"])
        for r in sorted(records, key=lambda x: x["id"]):
            w.writerow([r["id"], r["class"], r["_lbin"], r["_rep"], assign[r["_rep"]]])

    # ── summary ─────────────────────────────────────────────────────────────
    print("[4/4] summary")
    with open(os.path.join(splits_dir, "cv_summary.txt"), "w", encoding="utf-8") as f:
        f.write("CV splits summary\n=================\n\n")
        f.write(f"Sequences: {total} (pos {tot_pos}, neg {tot_neg})\n\n")
        f.write(f"[Homology-aware] {args.folds} folds, {len(clusters)} clusters\n")
        for fo in folds:
            f.write(f"  {fo['name']}: total={fo['n']} pos={fo['pos']} neg={fo['neg']}\n")
        f.write(f"\n[LOPO genus] {n_lopo_groups} folds (min_group={args.min_group})\n")
        gc = Counter(group.values())
        gpos = Counter(group[r["id"]] for r in records if r["class"] == "positive")
        for grp, c in gc.most_common():
            f.write(f"  {grp}: total={c} pos={gpos.get(grp,0)}\n")

    print(f"  cv_homology.tsv  ({args.folds} folds)")
    print(f"  cv_lopo.tsv      ({n_lopo_groups} folds)")
    for fo in folds:
        print(f"    {fo['name']}: total={fo['n']} pos={fo['pos']} neg={fo['neg']}")
    print(f"  summary: {os.path.join(splits_dir, 'cv_summary.txt')}")


if __name__ == "__main__":
    main()
