#!/usr/bin/env python3
# Matches each positive with negative(s) drawn from the real pool built in step
# 04, by TAXON -> LINEAGE -> LENGTH. Two details that matter:
#   - the negatives' taxid is read from 'OX=' in the UniProt header
#   - the positive's taxid is REMAPPED to 'proteome_taxid' (proteome_plan.tsv,
#     step 03), so IEDB internal ids (>=10^7) and resolved taxa can still match
#     by same_taxid against their real proteome
#
# Matching priority, searched over length bins [bin-1, bin, bin+1]:
#   1. same_taxid   (same effective taxon)
#   2. same_lineage (most specific ancestor that has negatives)
#   3. length_only  (anything of that length)
# Sampling is without replacement. R negatives per positive (--ratio).
#
# Needs ete3. No network.
#
#   python 05_match_negatives.py --species-dir species/Gallus_gallus --ratio 1
#
# in : <species-dir>/input/metadata.tsv
#      <species-dir>/negative_pool/candidates.no_homology.fasta
#      <species-dir>/negative_pool/proteome_plan.tsv   (taxid remap)
# out: <species-dir>/negative_pool/matched_negatives.fasta
#      <species-dir>/negative_pool/matching.tsv
#      <species-dir>/curated/sequences.tsv
#      <species-dir>/curated/fasta/all.fasta

import argparse
import csv
import hashlib
import os
import random
import re
from collections import defaultdict

from ete3 import NCBITaxa


def read_fasta(path):
    seqs, headers = {}, {}
    sid = header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(chunks)
                    headers[sid] = header
                header = line[1:]
                sid = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if sid is not None:
        seqs[sid] = "".join(chunks)
        headers[sid] = header
    return seqs, headers


def write_fasta(rows, out_file):
    with open(out_file, "w") as f:
        for r in rows:
            f.write(f">{r['id']}\n")
            s = r["sequence"]
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(rows, fields, out_file):
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in fields})


def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def first_taxid(raw):
    if not raw:
        return None
    for t in str(raw).replace("|", " ").split():
        t = t.strip()
        if t and t != "NA":
            return to_int(t)
    return None


def seq_md5(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def length_bin_label(n):
    if n is None:
        return "L00_missing"
    if n <= 15:  return "L01_very_short"
    if n <= 30:  return "L02_short_peptide"
    if n <= 60:  return "L03_long_peptide"
    if n <= 150: return "L04_short_protein"
    if n <= 400: return "L05_medium_protein"
    return "L06_long_protein"


def get_ox_taxid(header):
    m = re.search(r"OX=(\d+)", header)
    if m:
        return int(m.group(1))
    m = re.search(r"TaxID=(\d+)", header)         # in case the pool is UniRef
    return int(m.group(1)) if m else None


def get_os_name(header):
    m = re.search(r"OS=(.*?)(?:\s+[A-Z][A-Z]=|$)", header)
    return m.group(1).strip() if m else "NA"


def length_bin_number(length, bin_size):
    return (length - 1) // bin_size


def nearby_bins(center, n):
    bins = [center]
    for d in range(1, n + 1):
        bins += [center - d, center + d]
    return [b for b in bins if b >= 0]


def get_lineage(ncbi, taxid, cache):
    if taxid is None:
        return []
    if taxid not in cache:
        try:
            cache[taxid] = ncbi.get_lineage(taxid)
        except Exception:  # noqa: BLE001
            cache[taxid] = []
    return cache[taxid]


def choose(ids, available):
    ids = [x for x in ids if x in available]
    if not ids:
        return None
    c = random.choice(ids)
    available.remove(c)
    return c


def load_taxid_remap(plan_path):
    """original taxid -> proteome_taxid (real NCBI), for the resolved ones."""
    remap = {}
    if not os.path.isfile(plan_path):
        return remap
    for row in read_tsv(plan_path):
        if row.get("status") in ("resolved", "uniprotkb_taxon"):
            pt = to_int(row.get("proteome_taxid"))
            ot = to_int(row.get("taxid"))
            if ot is not None and pt is not None:
                remap[ot] = pt
    return remap


def main():
    ap = argparse.ArgumentParser(description="Match positives with negatives from the real pool")
    ap.add_argument("--species-dir", required=True)
    ap.add_argument("--ratio", type=int, default=1, help="negatives per positive")
    ap.add_argument("--bin-size", type=int, default=15)
    ap.add_argument("--neighbor-bins", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()
    random.seed(args.seed)

    sd = args.species_dir
    meta_file = os.path.join(sd, "input", "metadata.tsv")
    pool_file = os.path.join(sd, "negative_pool", "candidates.no_homology.fasta")
    plan_file = os.path.join(sd, "negative_pool", "proteome_plan.tsv")

    out_neg_fasta = os.path.join(sd, "negative_pool", "matched_negatives.fasta")
    out_matching = os.path.join(sd, "negative_pool", "matching.tsv")
    curated_dir = os.path.join(sd, "curated")
    fasta_dir = os.path.join(curated_dir, "fasta")
    out_sequences = os.path.join(curated_dir, "sequences.tsv")
    out_all_fasta = os.path.join(fasta_dir, "all.fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    for f in (meta_file, pool_file):
        if not os.path.isfile(f):
            raise SystemExit(f"Not found: {f}")

    metadata = read_tsv(meta_file)
    fields = list(metadata[0].keys())
    if "source" not in fields:
        fields.append("source")

    remap = load_taxid_remap(plan_file)
    print(f"[1/5] metadata: {len(metadata)} rows | taxid remap: {len(remap)} entries")

    # positives
    positives = []
    for row in metadata:
        if row.get("class") != "positive":
            continue
        if row.get("sequence_available") == "FALSE":
            continue
        seq = row.get("sequence", "")
        if not seq or seq == "NA":
            continue
        length = to_int(row.get("length")) or len(seq)
        ot = first_taxid(row.get("taxid"))
        eff = remap.get(ot, ot)                       # remapped to proteome_taxid
        positives.append({"id": row["id"], "taxid": eff, "length": length})
    print(f"[2/5] usable positives: {len(positives)}")

    # negative pool
    neg_pool, neg_headers = read_fasta(pool_file)
    ncbi = NCBITaxa()
    cache = {}
    neg_info = {}
    by_taxid_bin = defaultdict(list)
    by_lineage_bin = defaultdict(list)
    by_bin = defaultdict(list)
    for nid, seq in neg_pool.items():
        taxid = get_ox_taxid(neg_headers[nid])
        if taxid is None:
            continue
        length = len(seq)
        b = length_bin_number(length, args.bin_size)
        lineage = get_lineage(ncbi, taxid, cache)
        neg_info[nid] = {"taxid": taxid, "length": length, "bin": b}
        by_taxid_bin[(taxid, b)].append(nid)
        by_bin[b].append(nid)
        for anc in lineage:
            by_lineage_bin[(anc, b)].append(nid)
    available = set(neg_info)
    print(f"[3/5] pool: {len(neg_pool)} | indexable (with OX=): {len(available)}")

    # match
    print(f"[4/5] matching (ratio 1:{args.ratio})")
    random.shuffle(positives)
    matching_rows = []
    selected = []
    for pos in positives:
        pid, ptax, plen = pos["id"], pos["taxid"], pos["length"]
        pbin = length_bin_number(plen, args.bin_size)
        bins = nearby_bins(pbin, args.neighbor_bins)
        for _ in range(args.ratio):
            chosen, level, matched_taxon = None, "failed", "NA"
            if ptax is not None:
                for b in bins:
                    chosen = choose(by_taxid_bin[(ptax, b)], available)
                    if chosen:
                        level, matched_taxon = "same_taxid", str(ptax)
                        break
            if chosen is None and ptax is not None:
                for anc in reversed(get_lineage(ncbi, ptax, cache)):
                    for b in bins:
                        chosen = choose(by_lineage_bin[(anc, b)], available)
                        if chosen:
                            level, matched_taxon = "same_lineage", str(anc)
                            break
                    if chosen:
                        break
            if chosen is None:
                for b in bins:
                    chosen = choose(by_bin[b], available)
                    if chosen:
                        level = "length_only"
                        break
            if chosen:
                selected.append(chosen)
                ntax = neg_info[chosen]["taxid"]
                nlen = neg_info[chosen]["length"]
                diff = abs(plen - nlen)
            else:
                ntax = nlen = diff = "NA"
            matching_rows.append({
                "positive_id": pid, "positive_taxid": ptax if ptax is not None else "NA",
                "positive_length": plen, "negative_id": chosen or "NA",
                "negative_taxid": ntax, "negative_length": nlen,
                "length_diff": diff, "match_level": level, "matched_taxon": matched_taxon,
            })

    selected = sorted(set(selected))
    n_fail = sum(r["match_level"] == "failed" for r in matching_rows)
    levels = defaultdict(int)
    for r in matching_rows:
        levels[r["match_level"]] += 1
    print(f"  negatives selected: {len(selected)} | failed: {n_fail}")
    print(f"  levels: {dict(levels)}")

    # write
    print("[5/5] writing results")
    write_fasta([{"id": n, "sequence": neg_pool[n]} for n in selected], out_neg_fasta)
    write_tsv(matching_rows,
              ["positive_id", "positive_taxid", "positive_length", "negative_id",
               "negative_taxid", "negative_length", "length_diff", "match_level", "matched_taxon"],
              out_matching)

    sequence_rows = []
    for row in metadata:
        r = dict(row)
        r["source"] = "input_metadata"
        sequence_rows.append(r)
    for nid in selected:
        header = neg_headers[nid]
        seq = neg_pool[nid]
        length = len(seq)
        r = {c: "NA" for c in fields}
        r.update({
            "id": nid, "class": "negative", "label": "neg_matched",
            "confidence_level": "probable", "type_group": "neg_matched",
            "organism": get_os_name(header), "taxid": str(get_ox_taxid(header) or "NA"),
            "host_simple": metadata[0].get("host_simple", "NA"),
            "sequence": seq, "sequence_md5": seq_md5(seq), "length": str(length),
            "length_bin": length_bin_label(length), "references": "UniProt_reference_proteome",
            "sequence_available": "TRUE", "source": "matched_negative",
        })
        sequence_rows.append(r)
    write_tsv(sequence_rows, fields, out_sequences)
    write_fasta([{"id": r["id"], "sequence": r["sequence"]}
                 for r in sequence_rows if r.get("sequence") not in ("", "NA", None)],
                out_all_fasta)

    print(f"\nOK | total sequences: {len(sequence_rows)} "
          f"(pos {sum(1 for r in sequence_rows if r['class']=='positive')}, "
          f"neg {len(selected)})")
    print(f"  {out_sequences}")
    print(f"  {out_all_fasta}")


if __name__ == "__main__":
    main()
