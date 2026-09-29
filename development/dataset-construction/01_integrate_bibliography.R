#!/usr/bin/env Rscript
# Merges a host's manually curated BIBLIOGRAPHIC table
# (sequences_bibliographic.tsv) into the IEDB positives mined by step 00.
#
#   - every bibliographic record enters as a POSITIVE; the 'intensity' column
#     is not used
#       confidence_level = "bibliographic"
#       label            = "POS_SOURCE_ANTIGEN_BIBLIO"
#       type_group       = "source_antigen"
#   - same length filter as step 00: >= MIN_LENGTH (40 aa)
#   - 'pathogen' (a name) is resolved to an NCBI TAXID through the UniProt
#     taxonomy and cached, so these entries can get organism-matched negatives
#   - dedup by sequence MD5 against IEDB; on a clash the IEDB row wins
#
# Needs INTERNET for the taxid resolution.
#
#   Rscript scripts/01_integrate_bibliography.R \
#     --species-dir species/Gallus_gallus \
#     --bib "Gallus gallus/sequences_bibliographic.tsv"
#
# in : <species-dir>/input/metadata.tsv (from step 00) + the --bib TSV
# out: <species-dir>/input/metadata.tsv and pos.fasta, both updated;
#      the previous versions are backed up to input/backup_pre_biblio/

suppressMessages({
  library(tidyverse)
  library(digest)
  library(httr)
  library(jsonlite)
})

# ── Arguments ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
SPECIES_DIR <- get_arg("--species-dir")
BIB_TSV     <- get_arg("--bib")
MIN_LENGTH  <- as.integer(get_arg("--min-length", "40"))
CACHE_FILE  <- get_arg("--cache", "pathogen_taxid_cache.tsv")

if (is.null(SPECIES_DIR) || is.null(BIB_TSV)) {
  stop("Usage: Rscript 01_integrate_bibliography.R --species-dir species/<host> --bib <sequences_bibliographic.tsv>")
}

IEDB_META <- file.path(SPECIES_DIR, "input", "metadata.tsv")
OUT_META  <- IEDB_META
OUT_POS   <- file.path(SPECIES_DIR, "input", "pos.fasta")
BACKUP    <- file.path(SPECIES_DIR, "input", "backup_pre_biblio")

for (f in c(IEDB_META, BIB_TSV)) if (!file.exists(f)) stop(sprintf("Not found: %s", f))

# ── Helpers ─────────────────────────────────────────────────────────────────
is_canonical_aa <- function(x) !is.na(x) & x != "" & str_detect(x, "^[ACDEFGHIKLMNPQRSTVWY]+$")
md5_seq <- function(seq) vapply(seq, digest, FUN.VALUE = character(1), algo = "md5")

length_bin <- function(len) case_when(
  is.na(len) ~ "L00_missing",
  len <= 15  ~ "L01_very_short",
  len <= 30  ~ "L02_short_peptide",
  len <= 60  ~ "L03_long_peptide",
  len <= 150 ~ "L04_short_protein",
  len <= 400 ~ "L05_medium_protein",
  TRUE       ~ "L06_long_protein"
)

write_fasta <- function(df, path) {
  lines <- df %>%
    filter(sequence_available %in% c(TRUE, "TRUE")) %>%
    transmute(header = paste0(
      ">", id, "|class=", class, "|label=", label,
      "|confidence=", confidence_level, "|type=", type_group,
      "|length=", length, "|host=\"", host_simple, "\"",
      "|organism=\"", organism, "\"",
      "|taxid=", if_else(is.na(taxid) | taxid == "", "NA", as.character(taxid))),
      seq = sequence) %>%
    pivot_longer(c(header, seq), values_to = "line") %>% pull(line)
  writeLines(lines, path)
  cat(path, "->", nrow(df %>% filter(sequence_available %in% c(TRUE, "TRUE"))), "sequences\n")
}

# ── pathogen -> taxid (UniProt taxonomy, cached) ────────────────────────────
.tax_cache <- new.env()
if (file.exists(CACHE_FILE)) {
  cc <- suppressWarnings(read_tsv(CACHE_FILE, show_col_types = FALSE,
                                  col_types = cols(.default = col_character())))
  for (i in seq_len(nrow(cc))) assign(cc$pathogen[i], cc$taxid[i], envir = .tax_cache)
}

taxonomy_search <- function(query) {
  url <- paste0("https://rest.uniprot.org/taxonomy/search?format=json&size=5&query=",
                URLencode(query, reserved = TRUE))
  resp <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  out <- tryCatch(fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                           simplifyDataFrame = TRUE)$results,
                  error = function(e) NULL)
  out
}

name_to_taxid <- function(pathogen) {
  if (is.na(pathogen) || pathogen %in% c("", "NA")) return(NA_character_)
  if (exists(pathogen, envir = .tax_cache, inherits = FALSE))
    return(get(pathogen, envir = .tax_cache, inherits = FALSE))

  base <- str_squish(str_replace(pathogen, "\\(.*", ""))   # cut at the first parenthesis
  cands <- unique(c(pathogen, base))
  toks <- str_split(base, "\\s+")[[1]]
  while (length(toks) > 2) { toks <- head(toks, -1); cands <- unique(c(cands, paste(toks, collapse = " "))) }

  taxid <- NA_character_
  for (q in cands) {
    res <- taxonomy_search(q); Sys.sleep(0.2)
    if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) next
    sn <- tolower(res$scientificName %||% rep("", nrow(res)))
    hit <- which(str_detect(sn, fixed(tolower(q))))
    if (length(hit) >= 1) { taxid <- as.character(res$taxonId[hit[1]]); break }
    if (!is.null(res$rank)) {
      sp <- which(res$rank == "species")
      if (length(sp) >= 1) { taxid <- as.character(res$taxonId[sp[1]]); break }
    }
  }
  assign(pathogen, taxid, envir = .tax_cache)
  taxid
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Read the IEDB side ──────────────────────────────────────────────────────
cat("Reading IEDB metadata:", IEDB_META, "\n")
meta_iedb <- read_tsv(IEDB_META, show_col_types = FALSE, col_types = cols(.default = col_character()))
fields <- names(meta_iedb)
host_simple <- meta_iedb$host_simple[!is.na(meta_iedb$host_simple)][1]
md5_iedb <- unique(meta_iedb$sequence_md5[!is.na(meta_iedb$sequence_md5)])
cat(sprintf("  IEDB positives: %d | host: %s\n", sum(meta_iedb$class == "positive"), host_simple))

# ── Read and convert the bibliography ───────────────────────────────────────
cat("Reading bibliography:", BIB_TSV, "\n")
bib <- read_tsv(BIB_TSV, show_col_types = FALSE, col_types = cols(.default = col_character()))

bib <- bib %>%
  mutate(sequence = str_to_upper(str_replace_all(sequence, "\\s+", "")),
         length = nchar(sequence)) %>%
  filter(is_canonical_aa(sequence), length >= MIN_LENGTH)
cat(sprintf("  Bibliography after filtering (canonical, >=%d aa): %d\n", MIN_LENGTH, nrow(bib)))

# resolve a taxid per distinct pathogen
paths <- unique(bib$pathogen)
cat(sprintf("  Resolving %d pathogens to taxid...\n", length(paths)))
ptax <- tibble(pathogen = paths, taxid = vapply(paths, name_to_taxid, character(1)))
for (i in seq_len(nrow(ptax)))
  cat(sprintf("    %-40s -> %s\n", ptax$pathogen[i], ifelse(is.na(ptax$taxid[i]), "NA", ptax$taxid[i])))

# save the cache
all_cached <- tibble(pathogen = ls(.tax_cache),
                     taxid = vapply(ls(.tax_cache), function(k) get(k, envir = .tax_cache), character(1)))
write_tsv(all_cached, CACHE_FILE)

bib <- bib %>% left_join(ptax, by = "pathogen")

# build rows in the metadata schema
ext <- bib %>%
  mutate(
    class = "positive",
    confidence_level = "bibliographic",
    label = "POS_SOURCE_ANTIGEN_BIBLIO",
    candidate_test = FALSE, candidate_test_pos = FALSE, candidate_test_neg = FALSE,
    rel_clean = "source antigen", type_group = "source_antigen",
    organism = replace_na(pathogen, "NA"),
    host_simple = host_simple,
    sequence_md5 = md5_seq(sequence),
    length_bin = length_bin(length),
    stratify_primary = paste(class, type_group, length_bin, sep = "__"),
    stratify_secondary = paste(class, type_group, sep = "__"),
    pct_pos = NA, pct_subjects = NA, n_total = NA, n_pos = NA, n_neg = NA,
    references = replace_na(ref_id, "NA"),
    sequence_available = "TRUE",
    acc = coalesce(na_if(genbank_id, "NA"), na_if(uniprot_id, "NA"), na_if(original_id, "NA"), "NA"),
    id = paste0("EXT_", str_replace_all(acc, "[^A-Za-z0-9]", "_"), "_", substr(sequence_md5, 1, 12))
  )

# every schema column present, all as text, exactly as in the IEDB table
for (col in fields) if (!col %in% names(ext)) ext[[col]] <- NA
ext <- ext %>% select(all_of(fields)) %>% mutate(across(everything(), as.character))

# ── Dedup ─────────────────────────────────────────────────────────────────────
ext <- ext %>% distinct(sequence_md5, .keep_all = TRUE)        # within the bibliography
n_before <- nrow(ext)
ext_new <- ext %>% filter(!sequence_md5 %in% md5_iedb)         # against IEDB
cat(sprintf("  Unique bibliography: %d | duplicated in IEDB: %d | new: %d\n",
            n_before, n_before - nrow(ext_new), nrow(ext_new)))

n_sin_taxid <- sum(is.na(ext_new$taxid) | ext_new$taxid %in% c("", "NA"))
cat(sprintf("  EXT_ without taxid (will match by length only): %d\n", n_sin_taxid))

# ── Backup, merge, write ────────────────────────────────────────────────────
dir.create(BACKUP, showWarnings = FALSE, recursive = TRUE)
invisible(file.copy(IEDB_META, file.path(BACKUP, "metadata.tsv"), overwrite = TRUE))
if (file.exists(OUT_POS)) invisible(file.copy(OUT_POS, file.path(BACKUP, "pos.fasta"), overwrite = TRUE))

meta_comb <- bind_rows(meta_iedb, ext_new)
write_tsv(meta_comb, OUT_META)
cat(sprintf("\nCombined metadata: %d rows (%d positives)\n",
            nrow(meta_comb), sum(meta_comb$class == "positive")))

pos_comb <- meta_comb %>% filter(class == "positive")
write_fasta(pos_comb, OUT_POS)

cat("\n--- SUMMARY ---\n")
cat(sprintf("IEDB positives    : %d\n", sum(meta_iedb$class == "positive")))
cat(sprintf("EXT_ added        : %d\n", nrow(ext_new)))
cat(sprintf("Total positives   : %d\n", sum(meta_comb$class == "positive")))
cat("Done.\n")
