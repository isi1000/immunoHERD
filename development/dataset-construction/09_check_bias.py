#!/usr/bin/env python3
# Bias audit of the curated dataset (positives vs negatives) and of the CV
# folds. Looks for shortcuts the model could learn instead of immunogenicity:
#
#   - length pos vs neg (KS + per bin)
#   - amino-acid composition pos vs neg (total variation distance)
#   - organism/taxon: match_level breakdown and taxon overlap
#   - duplicate sequences and pos/neg conflicts (same MD5 in both classes)
#   - sequence leakage between homology folds; class balance per fold/group
#   - source of the positives (IEDB vs EXT_ bibliography)
#
# No network, read-only, stdlib only.
#
#   python 09_check_bias.py --species-dir species/Gallus_gallus
#
# out: the report, printed and saved to <sd>/qc/bias_report.txt

import argparse
import bisect
import csv
import os
from collections import Counter, defaultdict

AAS = "ACDEFGHIKLMNPQRSTVWY"


def read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def quantiles(vals, qs=(0, 0.25, 0.5, 0.75, 1.0)):
    if not vals:
        return {q: None for q in qs}
    s = sorted(vals)
    out = {}
    for q in qs:
        idx = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
        out[q] = s[idx]
    return out


def ks_stat(a, b):
    if not a or not b:
        return None
    a, b = sorted(a), sorted(b)
    na, nb = len(a), len(b)
    d = 0.0
    for v in sorted(set(a) | set(b)):
        d = max(d, abs(bisect.bisect_right(a, v) / na - bisect.bisect_right(b, v) / nb))
    return d


def aa_comp(seqs):
    """Mean per-sequence frequency of each AA (not weighted by length)."""
    acc = Counter()
    n = 0
    for s in seqs:
        L = len(s)
        if L == 0:
            continue
        c = Counter(s)
        for aa in AAS:
            acc[aa] += c.get(aa, 0) / L
        n += 1
    return ({aa: acc[aa] / n for aa in AAS} if n else {aa: 0.0 for aa in AAS}), n


def main():
    ap = argparse.ArgumentParser(description="Bias audit of the dataset")
    ap.add_argument("--species-dir", required=True)
    args = ap.parse_args()
    sd = args.species_dir

    seqs = read_tsv(os.path.join(sd, "curated", "sequences.tsv"))
    homo = read_tsv(os.path.join(sd, "splits", "cv_homology.tsv"))
    lopo = read_tsv(os.path.join(sd, "splits", "cv_lopo.tsv"))
    matching = read_tsv(os.path.join(sd, "negative_pool", "matching.tsv"))
    if not seqs:
        raise SystemExit("No curated/sequences.tsv")

    out = []
    def p(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    pos = [r for r in seqs if r.get("class") == "positive"]
    neg = [r for r in seqs if r.get("class") == "negative"]
    pos_seqs = [r.get("sequence", "") for r in pos if r.get("sequence") not in ("", "NA", None)]
    neg_seqs = [r.get("sequence", "") for r in neg if r.get("sequence") not in ("", "NA", None)]

    p("=" * 60)
    p(f"BIAS AUDIT — {sd}")
    p("=" * 60)

    # A. Class balance
    p("\n[A] Class balance")
    p(f"  positives: {len(pos)} | negatives: {len(neg)} | ratio neg/pos: {len(neg)/max(len(pos),1):.2f}")

    # B. Length
    p("\n[B] Length (aa) — pos vs neg")
    pl = [to_int(r.get("length")) or len(r.get("sequence", "")) for r in pos]
    nl = [to_int(r.get("length")) or len(r.get("sequence", "")) for r in neg]
    qp, qn = quantiles(pl), quantiles(nl)
    p(f"  pos  min/Q1/med/Q3/max: {qp[0]}/{qp[0.25]}/{qp[0.5]}/{qp[0.75]}/{qp[1.0]}")
    p(f"  neg  min/Q1/med/Q3/max: {qn[0]}/{qn[0.25]}/{qn[0.5]}/{qn[0.75]}/{qn[1.0]}")
    ks = ks_stat(pl, nl)
    flag = "  <-- CHECK (different distributions)" if ks and ks > 0.2 else "  OK"
    p(f"  KS(length) = {ks:.3f}{flag}" if ks is not None else "  KS: NA")
    p("  per length_bin (pos | neg):")
    bins = sorted(set(r.get("length_bin", "NA") for r in seqs))
    cb_pos = Counter(r.get("length_bin", "NA") for r in pos)
    cb_neg = Counter(r.get("length_bin", "NA") for r in neg)
    for b in bins:
        p(f"    {b:22s} {cb_pos.get(b,0):6d} | {cb_neg.get(b,0):6d}")

    # C. Amino-acid composition
    p("\n[C] Amino-acid composition — pos vs neg")
    cp, _ = aa_comp(pos_seqs)
    cn, _ = aa_comp(neg_seqs)
    tv = 0.5 * sum(abs(cp[a] - cn[a]) for a in AAS)
    p(f"  total variation distance (TV) = {tv:.4f} " +
      ("<-- CHECK (>0.05)" if tv > 0.05 else "OK"))
    diffs = sorted(((cp[a] - cn[a], a) for a in AAS), key=lambda x: abs(x[0]), reverse=True)
    p("  largest differences (AA: pos% vs neg%):")
    for d, a in diffs[:6]:
        p(f"    {a}: {cp[a]*100:5.2f}  vs {cn[a]*100:5.2f}   (Δ {d*100:+.2f})")

    # D. Organism / taxon (match_level + overlap)
    p("\n[D] Organism / taxon")
    if matching:
        lv = Counter(r.get("match_level", "NA") for r in matching)
        tot = sum(lv.values())
        for k in ("same_taxid", "same_lineage", "length_only", "failed"):
            p(f"  match_level {k:12s}: {lv.get(k,0):6d} ({100*lv.get(k,0)/max(tot,1):.1f}%)")
        lo = lv.get("length_only", 0) + lv.get("failed", 0)
        p("  -> " + ("CHECK: many negatives do not come from the positive's organism"
                     if lo / max(tot, 1) > 0.2 else "OK: most matched by taxon/lineage"))
    else:
        p("  (no matching.tsv)")
    pos_tax = set(r.get("taxid") for r in pos if r.get("taxid") not in ("", "NA", None))
    neg_tax = set(r.get("taxid") for r in neg if r.get("taxid") not in ("", "NA", None))
    p(f"  distinct taxa: pos={len(pos_tax)} neg={len(neg_tax)} shared={len(pos_tax & neg_tax)}")

    # E. Duplicates and pos/neg conflicts
    p("\n[E] Duplicate sequences and conflicts")
    md5_class = defaultdict(set)
    md5_count = Counter()
    for r in seqs:
        m = r.get("sequence_md5")
        if m and m != "NA":
            md5_class[m].add(r.get("class"))
            md5_count[m] += 1
    dups = sum(1 for m, c in md5_count.items() if c > 1)
    conflicts = [m for m, cls in md5_class.items() if "positive" in cls and "negative" in cls]
    p(f"  duplicate sequences by MD5 (>1 row): {dups}")
    p(f"  pos/neg CONFLICTS (same sequence in both classes): {len(conflicts)} " +
      ("<-- CRITICAL" if conflicts else "OK"))
    for m in conflicts[:5]:
        ids = [r["id"] for r in seqs if r.get("sequence_md5") == m]
        p(f"    md5 {m[:10]}…: {ids[:6]}")

    # F. Source of the positives
    p("\n[F] Source of the positives")
    src = Counter(("EXT_bibliography" if r["id"].startswith("EXT_") else "IEDB") for r in pos)
    for k, v in src.most_common():
        p(f"  {k}: {v}")

    # G. CV homology: balance and leakage
    p("\n[G] CV homology-aware")
    if homo:
        fold_of = {r["id"]: r.get("fold") for r in homo}
        cls_of = {r["id"]: r.get("class") for r in homo}
        per_fold = defaultdict(lambda: Counter())
        for r in homo:
            per_fold[r.get("fold")][r.get("class")] += 1
        for fo in sorted(per_fold):
            c = per_fold[fo]
            p(f"  {fo}: pos={c.get('positive',0)} neg={c.get('negative',0)}")
        # leakage: the same sequence (md5) in >1 fold
        md5_folds = defaultdict(set)
        id2md5 = {r["id"]: r.get("sequence_md5") for r in seqs}
        for r in homo:
            m = id2md5.get(r["id"])
            if m and m != "NA":
                md5_folds[m].add(r.get("fold"))
        leak = sum(1 for m, fs in md5_folds.items() if len(fs) > 1)
        p(f"  leakage (MD5 in >1 fold): {leak} " + ("<-- CHECK" if leak else "OK"))
    else:
        p("  (no cv_homology.tsv)")

    # H. CV LOPO: balance per group, and each genus in a single group
    p("\n[H] CV leave-one-pathogen-out (genus)")
    if lopo:
        per_grp = defaultdict(lambda: Counter())
        grp_genus = defaultdict(set)
        for r in lopo:
            per_grp[r.get("group")][r.get("class")] += 1
            grp_genus[r.get("group")].add(r.get("genus"))
        p(f"  groups (folds): {len(per_grp)}")
        for g in sorted(per_grp, key=lambda x: -sum(per_grp[x].values())):
            c = per_grp[g]
            p(f"    {g:28s} pos={c.get('positive',0):4d} neg={c.get('negative',0):5d}")
        # every genus should sit in exactly one group
        genus_groups = defaultdict(set)
        for r in lopo:
            if r.get("genus") not in ("NA", "", None):
                genus_groups[r.get("genus")].add(r.get("group"))
        split_genus = {g: grps for g, grps in genus_groups.items() if len(grps) > 1}
        p(f"  genera in >1 group (should be 0): {len(split_genus)} " +
          ("<-- CHECK" if split_genus else "OK"))
    else:
        p("  (no cv_lopo.tsv)")

    qc_dir = os.path.join(sd, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    rep = os.path.join(qc_dir, "bias_report.txt")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    p("\n" + "=" * 60)
    p(f"Report saved to: {rep}")


if __name__ == "__main__":
    main()
