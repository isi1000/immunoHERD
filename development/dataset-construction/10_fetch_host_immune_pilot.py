#!/usr/bin/env python3
# Pilot. Downloads from UniProt one marker immune gene per HOST (not per
# pathogen); B2M by default: beta-2-microglobulin, a universal MHC-I
# component, short and well annotated even in non-model species.
#
# The point is one sequence per host, embedded with the same ESM2 model as the
# epitopes, to get a real per-species "immune embedding" instead of a
# categorical host_id learnt from scratch out of six classes.
#
# Selection rules when a host has several candidate entries:
#   1. drop "(Fragment)" annotations (incomplete sequences)
#   2. prefer stronger protein evidence (lower protein_existence: 1 = protein
#      level, 5 = predicted)
#   3. on ties, keep the longest sequence (less risk of upstream truncation)
#   4. try reviewed:true (Swiss-Prot) first; if nothing, repeat without that
#      filter (TrEMBL)
#
# Needs INTERNET. stdlib urllib/json only, as in 03_download_proteomes.py.
#
#   python 10_fetch_host_immune_pilot.py --gene B2M --out-dir host_immune
#
# out: <out-dir>/<gene>_pilot.fasta          one sequence per host
#      <out-dir>/<gene>_pilot_provenance.tsv host, taxid, accession, length,
#                                            PE, reviewed

import argparse
import csv
import json
import os
import time
import urllib.parse
import urllib.request

UNIPROT = "https://rest.uniprot.org"
USER_AGENT = "immuno-model-thesis/1.0 (host immune pilot fetcher)"

# Same 6 hosts, same names, as everywhere else in the pipeline.
HOSTS_TAXA = {
    "Homo_sapiens":   9606,
    "Bos_taurus":     9913,
    "Sus_scrofa":     9823,
    "Equus_caballus": 9796,
    "Canis_sp":       9615,   # Canis lupus familiaris
    "Gallus_gallus":  9031,
}


def _fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


# Fallback for when the gene-symbol search returns nothing: some entries,
# mostly in non-model species, have the right ortholog annotated by protein
# NAME with the gene symbol field left empty. TAP2 in Equus_caballus is one:
# F7DH87 is the correct ortholog but carries no "gene". Manual map, only for
# genes already checked by hand.
PROTEIN_NAME_FALLBACK = {
    "B2M": "beta-2-microglobulin",
    "TAP1": "antigen peptide transporter 1",
    "TAP2": "antigen peptide transporter 2",
}


def _search_uniprot(gene, taxid, reviewed):
    query = f"gene:{gene} AND organism_id:{taxid}"
    if reviewed:
        query += " AND reviewed:true"
    return _run_search(query)


def _search_uniprot_by_name(protein_name, taxid):
    query = f'protein_name:"{protein_name}" AND organism_id:{taxid}'
    return _run_search(query)


def _run_search(query):
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,sequence,protein_existence,protein_name",
        "size": 20,
    }
    url = f"{UNIPROT}/uniprotkb/search?{urllib.parse.urlencode(params)}"
    data = _fetch_json(url)
    return data.get("results", [])


NON_STANDARD_AA = set("BJOUXZ")  # UniProt ambiguity / non-standard codes


def _pe_rank(protein_existence):
    """'1: Evidence at protein level' -> 1 (UniProt already puts it in the string)."""
    try:
        return int(protein_existence.split(":", 1)[0])
    except (ValueError, AttributeError, IndexError):
        return 9


def _pick_best(results):
    """Drops fragments (flag='Fragment', not just the name text) and sequences
    with ambiguous / non-standard residues; of the rest, picks the best entry
    (lowest PE, then longest)."""
    candidates = []
    for r in results:
        if r.get("proteinDescription", {}).get("flag") == "Fragment":
            continue
        seq = r.get("sequence", {}).get("value", "")
        if not seq or (set(seq) & NON_STANDARD_AA):
            continue
        pe_rank = _pe_rank(r.get("proteinExistence", ""))
        candidates.append((pe_rank, -len(seq), r["primaryAccession"], seq, r.get("proteinExistence", "NA")))
    if not candidates:
        return None
    candidates.sort()
    pe_rank, neg_len, acc, seq, pe_label = candidates[0]
    return {"accession": acc, "sequence": seq, "length": -neg_len, "protein_existence": pe_label}


def fetch_one(gene, host, taxid):
    for reviewed in (True, False):
        results = _search_uniprot(gene, taxid, reviewed)
        best = _pick_best(results)
        if best is not None:
            best["reviewed"] = reviewed
            best["via"] = "gene_symbol"
            return best
        time.sleep(0.3)

    protein_name = PROTEIN_NAME_FALLBACK.get(gene.upper())
    if protein_name:
        results = _search_uniprot_by_name(protein_name, taxid)
        best = _pick_best(results)
        if best is not None:
            best["reviewed"] = False
            best["via"] = "protein_name_fallback"
            return best
    return None


def main():
    ap = argparse.ArgumentParser(description="Pilot: one marker immune gene per host, from UniProt")
    ap.add_argument("--gene", default="B2M", help="Gene symbol to search for (default: B2M)")
    ap.add_argument("--out-dir", default="host_immune", help="Output folder")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    fasta_path = os.path.join(args.out_dir, f"{args.gene.lower()}_pilot.fasta")
    prov_path = os.path.join(args.out_dir, f"{args.gene.lower()}_pilot_provenance.tsv")

    rows = []
    with open(fasta_path, "w", encoding="utf-8") as fasta_out:
        for host, taxid in HOSTS_TAXA.items():
            best = fetch_one(args.gene, host, taxid)
            if best is None:
                print(f"[MISSING] {host} (taxid {taxid}): no usable entry for {args.gene}")
                rows.append({"host": host, "taxid": taxid, "accession": "NA", "length": "NA",
                             "protein_existence": "NA", "reviewed": "NA", "via": "NA"})
                continue
            fasta_out.write(f">{host}  ({args.gene}, UniProt {best['accession']})\n")
            fasta_out.write(best["sequence"] + "\n")   # no wrapping: these are short
            rows.append({"host": host, "taxid": taxid, "accession": best["accession"],
                         "length": best["length"], "protein_existence": best["protein_existence"],
                         "reviewed": best["reviewed"], "via": best["via"]})
            flag = "" if best["via"] == "gene_symbol" else "  [protein-name fallback]"
            print(f"[OK] {host:16s} {best['accession']:12s} len={best['length']:4d}  "
                  f"PE={best['protein_existence']}  reviewed={best['reviewed']}{flag}")

    with open(prov_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["host", "taxid", "accession", "length",
                                          "protein_existence", "reviewed", "via"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\n[OK] FASTA:      {fasta_path}")
    print(f"[OK] Provenance: {prov_path}")


if __name__ == "__main__":
    main()
