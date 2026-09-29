#!/usr/bin/env python3
# Resolves every antigen taxon with >= K positives to its UniProt REFERENCE
# PROTEOME and downloads the FASTA into a SHARED cache.
#
#   - strain -> species/genus resolution with ete3: a strain rarely has a
#     reference proteome of its own, its species does
#   - cache shared by UPID, so a proteome already there is never re-downloaded
#   - resolution cache (taxid -> UPID), so API calls are not repeated across
#     hosts and reruns
#
# Needs INTERNET, so it has to run on a node with network. ete3 for lineage,
# stdlib urllib/json/gzip for the rest; no 'requests'.
#
# The long tail (taxa with < K positives) is NOT downloaded here: step 05
# matches it by lineage/length against the proteomes already fetched.
#
#   python 03_download_proteomes.py --species-dir species/Gallus_gallus --min-pos 3
#
# in : <species-dir>/negative_pool/source_taxa.tsv   (from step 02)
# out: databases/proteomes/<UPID>.fasta.gz           proteome FASTA (shared)
#      databases/proteomes/_resolution.tsv           taxid -> UPID (shared)
#      <species-dir>/negative_pool/proteome_plan.tsv     per-taxon result
#      <species-dir>/negative_pool/proteomes_needed.txt  what this host needs

import argparse
import csv
import gzip
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

from ete3 import NCBITaxa

UNIPROT = "https://rest.uniprot.org"
USER_AGENT = "immuno-model-thesis/1.0 (proteome downloader)"
REFERENCE_TYPE = "Reference proteome"           # proteomeType value in the JSON
REFERENCE_FILTER = "proteome_type:REFERENCE"    # query facet, uppercase


# ── HTTP ────────────────────────────────────────────────────────────────────

def http_get(url, retries=4, timeout=120):
    """GET with retries. Returns bytes, or raises the last exception."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            wait = 2 ** attempt
            print(f"    warning: GET failed ({e}); retrying in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise last


def reference_proteomes_for_taxid(taxid):
    """Reference proteomes under an organism_id (subtree), filtered server-side.

    Uses the bounded 'search' endpoint plus proteome_type:REFERENCE, so a large
    taxon (E. coli and friends) does not drag thousands of proteomes back.
    """
    params = urllib.parse.urlencode({
        "format": "json",
        "size": "10",
        # taxonomy_id searches the subtree, so the reference is found even when
        # it is registered under a strain taxid below the one queried.
        "query": f"(taxonomy_id:{taxid}) AND ({REFERENCE_FILTER})",
    })
    url = f"{UNIPROT}/proteomes/search?{params}"
    data = http_get(url)
    try:
        return json.loads(data).get("results", [])
    except json.JSONDecodeError:
        return []


def pick_reference(results):
    """Picks the reference proteome; on ties, the largest protein_count."""
    if not results:
        return None
    refs = [r for r in results if r.get("proteomeType") == REFERENCE_TYPE] or results
    return sorted(refs, key=lambda r: r.get("proteinCount", 0) or 0, reverse=True)[0]


# ── Fallback by NAME, for IEDB internal taxids (>= 10^7) ────────────────────

def taxonomy_search(name):
    """Searches an organism name in the UniProt taxonomy."""
    params = urllib.parse.urlencode({"format": "json", "size": "5", "query": name})
    url = f"{UNIPROT}/taxonomy/search?{params}"
    try:
        return json.loads(http_get(url)).get("results", [])
    except Exception:  # noqa: BLE001
        return []


def name_to_taxid(organism):
    """Resolves an organism name to a real NCBI taxid.

    IEDB names usually carry the strain or isolate ('Newcastle disease virus
    (strain La Sota)', 'Pasteurella multocida X-73'), so progressively more
    general names are tried, and only a result whose scientificName CONTAINS
    the query is accepted, to avoid landing on the wrong clade.
    """
    organism = (organism or "").split("|")[0].strip()   # first organism if several
    if not organism:
        return None

    base = organism.split("(")[0].strip()               # cut at the first parenthesis
    cands = []
    for c in (organism, base):
        if c and c not in cands:
            cands.append(c)
    toks = base.split()                                 # drop strain/isolate suffixes
    while len(toks) > 2:
        toks = toks[:-1]
        c = " ".join(toks)
        if c not in cands:
            cands.append(c)

    for q in cands:
        results = taxonomy_search(q)
        time.sleep(0.2)
        ql = q.lower()
        for r in results:                               # 1) name containing the query
            sn = (r.get("scientificName") or "").lower()
            if ql and ql in sn:
                return str(r.get("taxonId"))
        for r in results:                               # 2) otherwise, anything at species rank
            if r.get("rank") == "species":
                return str(r.get("taxonId"))
    return None


# ── Taxon -> reference proteome ─────────────────────────────────────────────

def candidate_taxids(ncbi, taxid):
    """The taxon itself, then its species ancestor, then its genus."""
    cands = [taxid]
    try:
        lin = ncbi.get_lineage(int(taxid))
    except Exception:  # noqa: BLE001
        return cands
    if not lin:
        return cands
    ranks = ncbi.get_rank(lin)
    for want in ("species", "genus"):
        for t in reversed(lin):           # most specific first
            if ranks.get(t) == want and str(t) != str(taxid):
                cands.append(str(t))
                break
    return cands


def resolve_proteome(ncbi, taxid, res_cache):
    """Returns dict(upid, proteome_taxid, protein_count, level, status).

    Only successful resolutions are cached; 'no_proteome' is retried on every
    run, since they are few and a transient failure should not be frozen in.
    """
    cached = res_cache.get(taxid)
    if cached and cached.get("status") == "resolved":
        return cached

    result = {"upid": "", "proteome_taxid": "", "protein_count": "",
              "level": "", "status": "no_proteome"}
    for level, cand in zip(("self", "species", "genus"), candidate_taxids(ncbi, taxid)):
        results = reference_proteomes_for_taxid(cand)
        time.sleep(0.2)
        ref = pick_reference(results)
        if ref:
            tax = ref.get("taxonomy", {}) or {}
            result = {
                "upid": ref.get("id", ""),
                "proteome_taxid": str(tax.get("taxonId", cand)),
                "protein_count": ref.get("proteinCount", "") or "",
                "level": level,
                "status": "resolved",
            }
            break

    if result["status"] == "resolved":
        res_cache[taxid] = result
    return result


def resolve_taxon(ncbi, taxid, organism, res_cache):
    """Resolves by taxid; on failure, retries by organism NAME."""
    cached = res_cache.get(taxid)
    if cached and cached.get("status") == "resolved":
        return cached

    res = resolve_proteome(ncbi, taxid, res_cache)
    if res["status"] == "resolved":
        return res

    nm_taxid = name_to_taxid(organism)
    if nm_taxid and nm_taxid != taxid:
        res2 = resolve_proteome(ncbi, nm_taxid, res_cache)
        if res2["status"] == "resolved":
            res = dict(res2)
            res["level"] = "name:" + res2["level"]
            res_cache[taxid] = res          # cache under the original taxid too
    return res


# ── Proteome download ───────────────────────────────────────────────────────

def download_proteome(upid, cache_dir):
    """Downloads the proteome FASTA to cache_dir/<upid>.fasta.gz if missing."""
    dest = os.path.join(cache_dir, f"{upid}.fasta.gz")
    if os.path.isfile(dest) and os.path.getsize(dest) > 0:
        return dest, "cached"
    q = urllib.parse.quote(f"(proteome:{upid})")
    url = f"{UNIPROT}/uniprotkb/stream?compressed=true&format=fasta&query={q}"
    data = http_get(url, timeout=600)
    tmp = dest + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    # check it is gzip and holds at least one sequence
    try:
        with gzip.open(tmp, "rt") as f:
            first = f.readline()
        if not first.startswith(">"):
            os.remove(tmp)
            return dest, "empty"
    except OSError:
        os.remove(tmp)
        return dest, "bad_gzip"
    os.replace(tmp, dest)
    return dest, "downloaded"


# ── UniProtKB fallback, for taxa with NO reference proteome ─────────────────

def ncbi_real_taxid(ncbi, taxid, organism):
    """Best real NCBI taxid for the UniProtKB fallback: its own if real,
    otherwise the one resolved from the organism name."""
    try:
        t = int(taxid)
        if t < 10_000_000 and ncbi.get_lineage(t):
            return t
    except Exception:  # noqa: BLE001
        pass
    nm = name_to_taxid(organism)
    try:
        return int(nm) if nm else None
    except (TypeError, ValueError):
        return None


def uniprotkb_count(taxid):
    """Number of UniProtKB proteins under a taxon (subtree)."""
    params = urllib.parse.urlencode({"query": f"(taxonomy_id:{taxid})", "size": "0", "format": "list"})
    url = f"{UNIPROT}/uniprotkb/search?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as r:
            total = r.headers.get("x-total-results")
        return int(total) if total else None
    except Exception:  # noqa: BLE001
        return None


def download_uniprotkb_by_taxon(taxid, cache_dir):
    """Downloads ALL UniProtKB proteins of a taxon (subtree) as the pool."""
    upid = f"TAX{taxid}"
    dest = os.path.join(cache_dir, f"{upid}.fasta.gz")
    if os.path.isfile(dest) and os.path.getsize(dest) > 0:
        return upid, dest, "cached"
    q = urllib.parse.quote(f"(taxonomy_id:{taxid})")
    url = f"{UNIPROT}/uniprotkb/stream?compressed=true&format=fasta&query={q}"
    data = http_get(url, timeout=600)
    tmp = dest + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    try:
        with gzip.open(tmp, "rt") as f:
            first = f.readline()
        if not first.startswith(">"):
            os.remove(tmp)
            return upid, dest, "empty"
    except OSError:
        os.remove(tmp)
        return upid, dest, "bad_gzip"
    os.replace(tmp, dest)
    return upid, dest, "downloaded"


# ── Resolution cache (TSV) ──────────────────────────────────────────────────

def load_res_cache(path):
    cache = {}
    if not os.path.isfile(path):
        return cache
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cache[row["taxid"]] = {
                "upid": row.get("upid", ""),
                "proteome_taxid": row.get("proteome_taxid", ""),
                "protein_count": row.get("protein_count", ""),
                "level": row.get("level", ""),
                "status": row.get("status", "no_proteome"),
            }
    return cache


def save_res_cache(path, cache):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxid", "upid", "proteome_taxid", "protein_count", "level", "status"])
        for taxid, r in cache.items():
            w.writerow([taxid, r["upid"], r["proteome_taxid"],
                        r["protein_count"], r["level"], r["status"]])


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Resolve and download reference proteomes for one host")
    ap.add_argument("--species-dir", required=True, help="species/<host> (with negative_pool/source_taxa.tsv)")
    ap.add_argument("--cache", default="databases/proteomes", help="shared proteome cache")
    ap.add_argument("--min-pos", type=int, default=3, help="K: download taxa with >= K positives")
    ap.add_argument("--no-uniprotkb-fallback", dest="ukb_fallback", action="store_false",
                    help="do not fall back to UniProtKB when a taxon has no reference proteome")
    ap.add_argument("--ukb-max", type=int, default=100000,
                    help="max proteins for the UniProtKB fallback (avoids huge downloads)")
    ap.set_defaults(ukb_fallback=True)
    args = ap.parse_args()

    taxa_path = os.path.join(args.species_dir, "negative_pool", "source_taxa.tsv")
    if not os.path.isfile(taxa_path):
        raise SystemExit(f"Not found: {taxa_path} (did you run step 02?)")

    cache_dir = args.cache
    os.makedirs(cache_dir, exist_ok=True)
    res_cache_path = os.path.join(cache_dir, "_resolution.tsv")
    res_cache = load_res_cache(res_cache_path)

    host_neg_dir = os.path.join(args.species_dir, "negative_pool")
    out_plan = os.path.join(host_neg_dir, "proteome_plan.tsv")
    out_needed = os.path.join(host_neg_dir, "proteomes_needed.txt")

    # read this host's taxa
    taxa_rows = []
    with open(taxa_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["n_pos"] = int(row["n_pos"])
            taxa_rows.append(row)

    primary = [r for r in taxa_rows if r["n_pos"] >= args.min_pos]
    print(f"Total taxa: {len(taxa_rows)}")
    print(f"Taxa to download (>= {args.min_pos} positives): {len(primary)}")
    print(f"Cache: {cache_dir}\n")

    ncbi = NCBITaxa()

    plan = []
    needed_upids = OrderedDict()
    n_dl = n_cached = n_noproteome = n_ukb = 0

    for i, r in enumerate(primary, 1):
        taxid = r["taxid"]
        res = resolve_taxon(ncbi, taxid, r.get("organisms", ""), res_cache)
        status = res["status"]
        upid = res["upid"]
        dl_status = ""

        if status == "resolved" and upid:
            dest, dl_status = download_proteome(upid, cache_dir)
            if dl_status == "downloaded":
                n_dl += 1
            elif dl_status == "cached":
                n_cached += 1
            needed_upids[upid] = True
        elif args.ukb_fallback:
            # no reference proteome -> pull the whole UniProtKB taxon
            real = ncbi_real_taxid(ncbi, taxid, r.get("organisms", ""))
            cnt = uniprotkb_count(real) if real else None
            if real and (cnt is None or 0 < cnt <= args.ukb_max):
                upid, dest, dl_status = download_uniprotkb_by_taxon(real, cache_dir)
                if dl_status in ("downloaded", "cached"):
                    n_ukb += 1 if dl_status == "downloaded" else 0
                    n_cached += 1 if dl_status == "cached" else 0
                    needed_upids[upid] = True
                    res = {"upid": upid, "proteome_taxid": str(real),
                           "protein_count": str(cnt or ""), "level": "uniprotkb",
                           "status": "uniprotkb_taxon"}
                    status = "uniprotkb_taxon"
                else:
                    n_noproteome += 1
            else:
                n_noproteome += 1
                if real and cnt and cnt > args.ukb_max:
                    print(f"    UniProtKB for {real} skipped: {cnt} > ukb_max={args.ukb_max}", file=sys.stderr)
        else:
            n_noproteome += 1

        plan.append({
            "taxid": taxid,
            "n_pos": r["n_pos"],
            "organism": r.get("organisms", ""),
            "upid": upid,
            "proteome_taxid": res["proteome_taxid"],
            "protein_count": res["protein_count"],
            "resolved_level": res["level"],
            "status": status,
            "download": dl_status,
        })
        print(f"  [{i}/{len(primary)}] taxid {taxid} (n={r['n_pos']}) -> "
              f"{upid or 'NA'} [{res['level'] or '-'}/{dl_status or status}]")

    # save the updated resolution cache
    save_res_cache(res_cache_path, res_cache)

    # this host's plan
    with open(out_plan, "w", newline="", encoding="utf-8") as f:
        cols = ["taxid", "n_pos", "organism", "upid", "proteome_taxid",
                "protein_count", "resolved_level", "status", "download"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(plan)

    # list of proteomes this host needs, for step 04
    with open(out_needed, "w", encoding="utf-8") as f:
        for upid in needed_upids:
            f.write(os.path.join(cache_dir, f"{upid}.fasta.gz") + "\n")

    print("\n--- Summary ---")
    print(f"Downloaded (reference) : {n_dl}")
    print(f"UniProtKB (fallback)   : {n_ukb}")
    print(f"Already cached         : {n_cached}")
    print(f"No proteome            : {n_noproteome}")
    print(f"Unique proteomes for this host: {len(needed_upids)}")
    print(f"\nPlan:     {out_plan}")
    print(f"Needed:   {out_needed}")
    print(f"Resolut.: {res_cache_path}")


if __name__ == "__main__":
    main()
