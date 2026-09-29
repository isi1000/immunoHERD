#!/usr/bin/env Rscript
# Mines an IEDB export for POSITIVES ONLY, per host, at protein/peptide level.
#
# It deliberately does three things:
#   1. no negatives of any kind are produced here; they are built later from
#      reference proteomes (steps 02-05)
#   2. only positives with a DIRECT relation, or with an assay organism, are
#      kept:
#         epitope                    -> POS_EPITOPE_STRONG / _WEAK
#         fragment of source antigen -> POS_FRAGMENT_STRONG / _WEAK
#         source antigen             -> POS_SOURCE_ANTIGEN_STRONG / _WEAK
#         source organism (is_pos)   -> POS_ASSAY_ANTIGEN_WEAK (always weak)
#      Indirect positives, rescued through a taxonomic parent/child or through
#      an empty relation, are dropped: too far from the measured antigen.
#   3. a minimum length of MIN_LENGTH (40 aa) applies to EVERY positive, which
#      removes almost all short epitopes and leaves a protein-level dataset
#
# Needs INTERNET (sequences are fetched from UniProt/NCBI by accession).
#
#   Rscript scripts/00_mine_iedb_positives.R
#
# out, per host, in OUTDIR_ROOT/<host>/input/:
#   metadata.tsv   positive metadata
#   pos.fasta      positive sequences
#   qc/*.tsv       QC: what was lost at each step
#   summary.tsv    host summary

library(tidyverse)
library(stringr)
library(httr)
library(digest)


# ============================================================
# 0. Configuration  <- EDIT HERE
# ============================================================

INPUT_CSV   <- "bcell_full_v3.csv"   # raw IEDB export
OUTDIR_ROOT <- "species_def"         # output root; writes to <root>/<host>/input/
MIN_ASSAYS  <- 400                   # minimum assays for a host to be processed
MIN_LENGTH  <- 40                    # minimum sequence length (aa)
SLEEP_SEC   <- 0.15                  # pause between sequence downloads

# Optional environment overrides, handy when rerunning one host:
#   OUTDIR_ROOT=...            change the output folder
#   IEDB_HOSTS="Host A,Host B" process only those hosts (still above MIN_ASSAYS)
if (nzchar(Sys.getenv("OUTDIR_ROOT"))) OUTDIR_ROOT <- Sys.getenv("OUTDIR_ROOT")


# ============================================================
# 1. Helpers
# ============================================================

is_canonical_aa <- function(x) {
  !is.na(x) & x != "" & str_detect(x, "^[ACDEFGHIKLMNPQRSTVWY]+$")
}

get_col <- function(pattern, names_df) {
  grep(pattern, names_df, value = TRUE)[1]
}

extract_accession_from_iri <- function(iri,
                                       parent_iri = NA_character_,
                                       epi_parent_iri = NA_character_) {
  case_when(
    str_detect(iri, "uniprot.org/uniprot/") ~
      str_extract(iri, "(?<=uniprot/)[A-Z0-9]+(?:-[0-9]+)?"),
    str_detect(iri, "ncbi.nlm.nih.gov/protein/") ~
      str_extract(iri, "(?<=protein/)[A-Z0-9_.]+"),
    !is.na(parent_iri) & str_detect(parent_iri, "uniprot.org/uniprot/") ~
      str_extract(parent_iri, "(?<=uniprot/)[A-Z0-9]+(?:-[0-9]+)?"),
    !is.na(parent_iri) & str_detect(parent_iri, "ncbi.nlm.nih.gov/protein/") ~
      str_extract(parent_iri, "(?<=protein/)[A-Z0-9_.]+"),
    !is.na(epi_parent_iri) & str_detect(epi_parent_iri, "uniprot.org/uniprot/") ~
      str_extract(epi_parent_iri, "(?<=uniprot/)[A-Z0-9]+(?:-[0-9]+)?"),
    !is.na(epi_parent_iri) & str_detect(epi_parent_iri, "ncbi.nlm.nih.gov/protein/") ~
      str_extract(epi_parent_iri, "(?<=protein/)[A-Z0-9_.]+"),
    TRUE ~ NA_character_
  )
}

is_uniprot_acc <- function(acc) {
  str_detect(
    acc,
    "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$"
  )
}

fetch_sequence <- function(acc) {
  if (is.na(acc) || acc == "") {
    return(tibble(accession = acc, sequence = NA_character_,
                  source = NA_character_, status = NA_integer_))
  }

  if (is_uniprot_acc(acc)) {
    url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".fasta")
    source <- "uniprot"
  } else {
    url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      "?db=protein&id=", acc, "&rettype=fasta&retmode=text"
    )
    source <- "ncbi"
  }

  resp <- tryCatch(GET(url), error = function(e) NULL)

  if (is.null(resp)) {
    return(tibble(accession = acc, sequence = NA_character_,
                  source = source, status = NA_integer_))
  }
  if (status_code(resp) != 200) {
    return(tibble(accession = acc, sequence = NA_character_,
                  source = source, status = status_code(resp)))
  }

  txt   <- content(resp, as = "text", encoding = "UTF-8")
  lines <- str_split(txt, "\n", simplify = TRUE) %>% as.character()
  lines <- lines[lines != ""]
  seq   <- paste0(lines[!str_starts(lines, ">")], collapse = "")

  tibble(accession = acc, sequence = seq, source = source, status = 200)
}

clean_folder_name <- function(x) {
  str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

length_bin <- function(len) {
  case_when(
    is.na(len) ~ "L00_missing",
    len <= 15  ~ "L01_very_short",
    len <= 30  ~ "L02_short_peptide",
    len <= 60  ~ "L03_long_peptide",
    len <= 150 ~ "L04_short_protein",
    len <= 400 ~ "L05_medium_protein",
    TRUE       ~ "L06_long_protein"
  )
}

write_fasta_final <- function(df, outfile) {
  fasta_lines <- df %>%
    filter(sequence_available) %>%
    transmute(
      header = paste0(
        ">", id,
        "|class=", class,
        "|label=", label,
        "|confidence=", confidence_level,
        "|type=", type_group,
        "|length=", length,
        "|host=\"", host_simple, "\"",
        "|organism=\"", organism, "\"",
        "|taxid=", if_else(is.na(taxid) | taxid == "", "NA", taxid)
      ),
      seq = sequence
    ) %>%
    pivot_longer(cols = c(header, seq), values_to = "line") %>%
    pull(line)

  writeLines(fasta_lines, outfile)
  cat(outfile, "->", length(fasta_lines) / 2, "sequences\n")
}


# ============================================================
# 2. Read the export and rebuild the column names
# ============================================================

raw <- read_delim(
  INPUT_CSV,
  delim = ",",
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

header <- raw[1, ]
df <- raw[-1, ]

names(df) <- map_chr(names(df), function(col) {
  h <- str_trim(as.character(header[[col]]))
  if (!is.na(h) && h != "") paste0(h, " (", col, ")") else col
})

cat("Total rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n\n")


# ============================================================
# 3. Locate the columns needed
# ============================================================

cols <- list(
  rel1                  = get_col("^Epitope Relation \\(1st immunogen", names(df)),
  obj1                  = get_col("^Object Type \\(1st immunogen", names(df)),
  name1                 = get_col("^Name \\(1st immunogen", names(df)),
  iri1                  = get_col("^IRI \\(1st immunogen", names(df)),
  org1                  = get_col("^Source Organism \\(1st immunogen", names(df)),
  org1_iri              = get_col("^Source Organism IRI \\(1st immunogen", names(df)),

  assay                 = get_col("^Qualitative Measure \\(Assay", names(df)),
  reference             = get_col("^IEDB IRI \\(Reference", names(df)),
  host                  = get_col("^Name \\(Host", names(df)),
  n_subjects_tested     = get_col("^Number of Subjects Tested \\(Assay", names(df)),
  n_subjects_positive   = get_col("^Number of Subjects Positive \\(Assay", names(df)),

  assay_source_mol_iri  = get_col("^Source Molecule IRI \\(Assay Antigen", names(df)),
  assay_parent_iri      = get_col("^Molecule Parent IRI \\(Assay Antigen", names(df)),
  assay_antigen_org     = get_col("^Source Organism \\(Assay Antigen", names(df)),
  assay_antigen_org_iri = get_col("^Source Organism IRI \\(Assay Antigen", names(df)),

  epi_parent_iri        = get_col("^Molecule Parent IRI \\(Epitope", names(df)),
  epi_start             = get_col("^Starting Position \\(Epitope", names(df)),
  epi_end               = get_col("^Ending Position \\(Epitope", names(df))
)

missing_cols <- names(cols)[map_lgl(cols, is.na)]
if (length(missing_cols) > 0) {
  cat("Columns not found:\n"); print(missing_cols)
  cat("\nAvailable columns:\n");  print(names(df))
  stop("Required columns are missing.")
}
cat("All required columns found.\n\n")


# ============================================================
# 4. Collapse host names and keep hosts with > MIN_ASSAYS assays
# ============================================================

df_hosts <- df %>%
  mutate(
    host_raw = str_squish(.data[[cols$host]]),
    host_simple = case_when(
      str_detect(host_raw, regex("Homo sapiens", ignore_case = TRUE)) ~ "Homo sapiens",
      str_detect(host_raw, regex("Mus musculus", ignore_case = TRUE)) ~ "Mus musculus",
      str_detect(host_raw, regex("Oryctolagus cuniculus", ignore_case = TRUE)) ~ "Oryctolagus cuniculus",
      str_detect(host_raw, regex("Rattus norvegicus|Rattus \\(rat\\)|Rattus rattus", ignore_case = TRUE)) ~ "Rattus sp.",
      str_detect(host_raw, regex("Macaca", ignore_case = TRUE)) ~ "Macaca sp.",
      str_detect(host_raw, regex("Aotus|Papio|Saimiri|Callithrix|Cebuella|Alouatta|Pithecia|Saguinus|Pan troglodytes|Chlorocebus|Platyrrhini|Simiiformes", ignore_case = TRUE)) ~ "Non-human primate",
      str_detect(host_raw, regex("Sus scrofa", ignore_case = TRUE)) ~ "Sus scrofa",
      str_detect(host_raw, regex("Bos taurus|Bos indicus|hybrid cow|bovine", ignore_case = TRUE)) ~ "Bos sp.",
      str_detect(host_raw, regex("Ovis aries", ignore_case = TRUE)) ~ "Ovis aries",
      str_detect(host_raw, regex("Capra hircus", ignore_case = TRUE)) ~ "Capra hircus",
      str_detect(host_raw, regex("Canis lupus|Canis familiaris|dogs|wolf or dog", ignore_case = TRUE)) ~ "Canis sp.",
      str_detect(host_raw, regex("Equus caballus|Equid", ignore_case = TRUE)) ~ "Equus caballus",
      str_detect(host_raw, regex("Gallus gallus", ignore_case = TRUE)) ~ "Gallus gallus",
      str_detect(host_raw, regex("Cavia porcellus", ignore_case = TRUE)) ~ "Cavia porcellus",
      str_detect(host_raw, regex("Lama glama|Vicugna pacos|Camelus|Camelid", ignore_case = TRUE)) ~ "Camelidae",
      str_detect(host_raw, regex("Mesocricetus auratus|Cricetinae|hamster", ignore_case = TRUE)) ~ "Hamster",
      str_detect(host_raw, regex("Mustela|mink", ignore_case = TRUE)) ~ "Mustelidae",
      str_detect(host_raw, regex("Anas|Cairina|duck|goose|Anser", ignore_case = TRUE)) ~ "Anseriformes",
      str_detect(host_raw, regex("Fish|salmon|tilapia|seabass|carp|halibut|trout|grouper|Dicentrarchus|Oreochromis|Oncorhynchus|Cyprinus|Paralichthys|Epinephelus|Morone|Lates", ignore_case = TRUE)) ~ "Fish",
      str_detect(host_raw, regex("shark|dogfish|Chiloscyllium|Squalus|Ginglymostoma|Orectolobus", ignore_case = TRUE)) ~ "Cartilaginous fish",
      TRUE ~ host_raw
    )
  )

host_counts <- df_hosts %>%
  count(host_simple, name = "n_assays") %>%
  arrange(desc(n_assays))

dir.create(OUTDIR_ROOT, showWarnings = FALSE, recursive = TRUE)
write.table(host_counts, file.path(OUTDIR_ROOT, "host_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hosts_keep <- host_counts %>% filter(n_assays > MIN_ASSAYS) %>% pull(host_simple)

# Optional host filter (IEDB_HOSTS="Host A,Host B")
host_filter <- Sys.getenv("IEDB_HOSTS", "")
if (nzchar(host_filter)) {
  wanted <- str_squish(str_split(host_filter, ",")[[1]])
  hosts_keep <- intersect(hosts_keep, wanted)
  cat("IEDB_HOSTS filter active ->", paste(hosts_keep, collapse = ", "), "\n")
}

cat("Hosts with more than", MIN_ASSAYS, "assays:\n")
print(host_counts %>% filter(n_assays > MIN_ASSAYS))
cat("\n")


# ============================================================
# 5. Process each host
# ============================================================

# DIRECT relations kept, at protein/peptide level.
# "source organism" is kept only when is_pos (POS_ASSAY_ANTIGEN_WEAK).
direct_relations <- c(
  "epitope",
  "fragment of source antigen",
  "source antigen",
  "source organism"
)

for (HOST in hosts_keep) {

  cat("\n============================================================\n")
  cat("Processing host:", HOST, "\n")
  cat("============================================================\n")

  outdir <- file.path(OUTDIR_ROOT, clean_folder_name(HOST), "input")
  qcdir  <- file.path(outdir, "qc")
  dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)

  df_host <- df_hosts %>% filter(host_simple == HOST)
  cat("Rows for", HOST, ":", nrow(df_host), "\n\n")


  # ----------------------------------------------------------
  # 5.1. Base table (direct relations only)
  # ----------------------------------------------------------

  df_base <- df_host %>%
    mutate(
      rel_clean = str_to_lower(str_squish(.data[[cols$rel1]])),
      obj_clean = str_to_lower(str_squish(.data[[cols$obj1]])),
      qm_clean  = str_to_lower(str_squish(.data[[cols$assay]])),

      is_pos = str_detect(qm_clean, "positive"),
      is_neg = str_detect(qm_clean, "negative"),

      n_subjects_tested   = as.numeric(.data[[cols$n_subjects_tested]]),
      n_subjects_positive = as.numeric(.data[[cols$n_subjects_positive]]),

      seq_1st = .data[[cols$name1]] %>%
        str_squish() %>% str_to_upper() %>% str_replace_all("\\s+", ""),

      acc_source_antigen = extract_accession_from_iri(
        iri = str_squish(.data[[cols$iri1]])
      ),
      acc_source_organism = extract_accession_from_iri(
        iri            = str_squish(.data[[cols$assay_source_mol_iri]]),
        parent_iri     = str_squish(.data[[cols$assay_parent_iri]]),
        epi_parent_iri = str_squish(.data[[cols$epi_parent_iri]])
      ),

      value = case_when(
        rel_clean == "epitope"                     ~ seq_1st,
        rel_clean == "fragment of source antigen"  ~ seq_1st,
        rel_clean == "source antigen"              ~ acc_source_antigen,
        rel_clean == "source organism" & is_pos    ~ acc_source_organism,
        TRUE ~ NA_character_
      ),
      value_type = case_when(
        rel_clean == "epitope"                     ~ "epitope_sequence",
        rel_clean == "fragment of source antigen"  ~ "fragment_sequence",
        rel_clean == "source antigen"              ~ "source_antigen_acc",
        rel_clean == "source organism" & is_pos    ~ "assay_antigen_acc",
        TRUE ~ "discard"
      ),

      organism = case_when(
        rel_clean == "source organism" & is_pos ~ str_squish(.data[[cols$assay_antigen_org]]),
        TRUE ~ str_squish(.data[[cols$org1]])
      ),
      organism_iri = case_when(
        rel_clean == "source organism" & is_pos ~ str_squish(.data[[cols$assay_antigen_org_iri]]),
        TRUE ~ str_squish(.data[[cols$org1_iri]])
      ),
      source_taxid = str_extract(
        organism_iri, "(?<=NCBITaxon_)[0-9]+|(?<=taxon/)[0-9]+"
      ),

      epi_start = as.numeric(.data[[cols$epi_start]]),
      epi_end   = as.numeric(.data[[cols$epi_end]]),
      reference = .data[[cols$reference]]
    ) %>%
    filter(rel_clean %in% direct_relations) %>%
    filter(obj_clean != "non-peptidic") %>%
    filter(value_type != "discard") %>%
    filter(!is.na(value), value != "")

  cat("Rows in df_base:", nrow(df_base), "\n\n")

  if (nrow(df_base) == 0) {
    cat("No usable data for", HOST, ". Skipping.\n")
    next
  }


  # ----------------------------------------------------------
  # 5.2. Summary per distinct immunogen
  # ----------------------------------------------------------

  immunogen_stats <- df_base %>%
    group_by(rel_clean, value_type, value) %>%
    summarise(
      n_total = n(),
      n_pos   = sum(is_pos, na.rm = TRUE),
      n_neg   = sum(is_neg, na.rm = TRUE),
      pct_pos = 100 * n_pos / n_total,

      sum_subjects_tested   = sum(n_subjects_tested, na.rm = TRUE),
      sum_subjects_positive = sum(n_subjects_positive, na.rm = TRUE),
      subjects_available    = any(!is.na(n_subjects_tested) & n_subjects_tested > 0),

      organism     = paste(sort(unique(na.omit(organism))), collapse = " | "),
      source_taxid = paste(sort(unique(na.omit(source_taxid))), collapse = " | "),
      references   = paste(sort(unique(na.omit(reference))), collapse = " | "),
      .groups = "drop"
    ) %>%
    mutate(
      pct_subjects = if_else(
        subjects_available & sum_subjects_tested > 0,
        100 * sum_subjects_positive / sum_subjects_tested,
        NA_real_
      )
    )


  # ----------------------------------------------------------
  # 5.3. Positive selection (strong / weak)
  # ----------------------------------------------------------
  # strong_positive: n_pos >= 2 and pct_pos >= 50, direct relations only
  #                  (epitope / fragment / source antigen)
  # weak_positive  : n_pos >= 1 and pct_pos >= 50, whatever is not strong.
  #                  "source organism" can only ever be weak.
  # ----------------------------------------------------------

  pos_strong <- immunogen_stats %>%
    filter(
      !is.na(rel_clean),
      rel_clean != "source organism",
      n_pos >= 2,
      pct_pos >= 50
    ) %>%
    mutate(
      class = "positive",
      confidence_level = "strong_positive",
      label = case_when(
        rel_clean == "epitope"                    ~ "POS_EPITOPE_STRONG",
        rel_clean == "fragment of source antigen" ~ "POS_FRAGMENT_STRONG",
        rel_clean == "source antigen"             ~ "POS_SOURCE_ANTIGEN_STRONG",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(label))

  pos_weak <- immunogen_stats %>%
    filter(n_pos >= 1, pct_pos >= 50) %>%
    anti_join(pos_strong, by = c("rel_clean", "value_type", "value")) %>%
    mutate(
      class = "positive",
      confidence_level = "weak_positive",
      label = case_when(
        rel_clean == "epitope"                    ~ "POS_EPITOPE_WEAK",
        rel_clean == "fragment of source antigen" ~ "POS_FRAGMENT_WEAK",
        rel_clean == "source antigen"             ~ "POS_SOURCE_ANTIGEN_WEAK",
        rel_clean == "source organism"            ~ "POS_ASSAY_ANTIGEN_WEAK",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(label))

  pos_selected <- bind_rows(pos_strong, pos_weak)

  if (nrow(pos_selected) == 0) {
    cat("No positives for", HOST, ". Skipping.\n")
    next
  }


  # ----------------------------------------------------------
  # 5.4. Fetch sequences by accession (source antigen / source organism)
  # ----------------------------------------------------------

  cat("--- Fetching sequences ---\n")

  all_unique_accessions <- pos_selected %>%
    filter(value_type %in% c("source_antigen_acc", "assay_antigen_acc")) %>%
    pull(value) %>% unique() %>% na.omit()

  cat("Unique accessions to fetch:", length(all_unique_accessions), "\n")

  if (length(all_unique_accessions) > 0) {
    global_sequence_catalog <- map_dfr(all_unique_accessions, function(acc) {
      Sys.sleep(SLEEP_SEC)
      fetch_sequence(acc)
    }) %>%
      filter(!is.na(sequence), sequence != "") %>%
      mutate(protein_length = nchar(sequence))
  } else {
    global_sequence_catalog <- tibble(
      accession = character(), sequence = character(),
      source = character(), status = integer(), protein_length = integer()
    )
  }

  cat("Sequences fetched:", nrow(global_sequence_catalog), "\n\n")


  # ----------------------------------------------------------
  # 5.5. Attach the final sequence to each positive
  # ----------------------------------------------------------
  # epitope / fragment               -> the sequence itself (value)
  # source antigen / source organism -> the sequence fetched by accession
  # ----------------------------------------------------------

  pos_selected <- pos_selected %>%
    left_join(
      global_sequence_catalog %>% select(accession, sequence),
      by = c("value" = "accession")
    ) %>%
    mutate(
      final_sequence = case_when(
        rel_clean == "epitope"                    ~ value,
        rel_clean == "fragment of source antigen" ~ value,
        rel_clean %in% c("source antigen", "source organism") ~ sequence,
        TRUE ~ NA_character_
      )
    ) %>%
    select(-sequence)


  # ----------------------------------------------------------
  # 5.6. Build the metadata
  # ----------------------------------------------------------

  selected_seq <- pos_selected %>%
    mutate(
      sequence_key = paste(HOST, rel_clean, value, class, label, sep = "|"),
      hash_id = substr(
        vapply(sequence_key, digest, FUN.VALUE = character(1), algo = "md5"), 1, 12
      ),
      id = case_when(
        rel_clean == "source antigen"             ~ paste0("SA_",  value, "_", hash_id),
        rel_clean == "source organism"            ~ paste0("AS_",  value, "_", hash_id),
        rel_clean == "epitope"                    ~ paste0("EPI_", hash_id),
        rel_clean == "fragment of source antigen" ~ paste0("FRAG_", hash_id),
        TRUE ~ paste0("UNK_", hash_id)
      )
    )

  metadata <- selected_seq %>%
    mutate(
      host_simple = HOST,
      sequence = str_to_upper(str_replace_all(final_sequence, "\\s+", "")),
      sequence_available = !is.na(sequence) & sequence != "",
      sequence_md5 = if_else(
        sequence_available,
        vapply(sequence, digest, FUN.VALUE = character(1), algo = "md5"),
        NA_character_
      ),
      length = nchar(sequence),

      type_group = case_when(
        rel_clean == "epitope"                    ~ "epitope",
        rel_clean == "fragment of source antigen" ~ "fragment",
        rel_clean == "source antigen"             ~ "source_antigen",
        rel_clean == "source organism"            ~ "assay_antigen",
        TRUE ~ "unknown"
      ),
      length_bin = length_bin(length),
      taxid = source_taxid,

      candidate_test     = confidence_level == "strong_positive",
      candidate_test_pos = confidence_level == "strong_positive",
      candidate_test_neg = FALSE,

      stratify_primary   = paste(class, type_group, length_bin, sep = "__"),
      stratify_secondary = paste(class, type_group, sep = "__")
    ) %>%
    select(
      id, class, label, confidence_level,
      candidate_test, candidate_test_pos, candidate_test_neg,
      rel_clean, type_group, organism, taxid, host_simple,
      sequence, sequence_md5, length, length_bin,
      stratify_primary, stratify_secondary,
      pct_pos, pct_subjects, n_total, n_pos, n_neg,
      references, sequence_available
    )


  # ----------------------------------------------------------
  # 5.7. Final filter: canonical, available, and length >= MIN_LENGTH
  # ----------------------------------------------------------

  metadata_before_filter <- metadata %>%
    mutate(
      canonical_sequence = is_canonical_aa(sequence),
      fail_no_sequence   = !sequence_available,
      fail_noncanonical  = !canonical_sequence,
      fail_too_short     = is.na(length) | length < MIN_LENGTH,
      pass_final_filter  = sequence_available & canonical_sequence &
                           !is.na(length) & length >= MIN_LENGTH,
      final_filter_reason = case_when(
        pass_final_filter ~ "PASS",
        fail_no_sequence  ~ "NO_SEQUENCE",
        fail_noncanonical ~ "NON_CANONICAL_SEQUENCE",
        fail_too_short    ~ "TOO_SHORT",
        TRUE              ~ "OTHER"
      )
    )

  metadata <- metadata_before_filter %>%
    filter(pass_final_filter) %>%
    select(-canonical_sequence, -fail_no_sequence, -fail_noncanonical,
           -fail_too_short, -pass_final_filter, -final_filter_reason)


  # ----------------------------------------------------------
  # 5.8. QC
  # ----------------------------------------------------------

  qc_selection_counts <- tibble(
    step = c(
      "01_df_host_rows",
      "02_df_base_rows",
      "03_immunogen_stats_unique",
      "04_pos_strong",
      "05_pos_weak",
      "06_pos_selected",
      "07_metadata_before_final_filter",
      "08_metadata_after_final_filter",
      "09_final_strong",
      "10_final_weak"
    ),
    n_rows = c(
      nrow(df_host),
      nrow(df_base),
      nrow(immunogen_stats),
      nrow(pos_strong),
      nrow(pos_weak),
      nrow(pos_selected),
      nrow(metadata_before_filter),
      nrow(metadata),
      sum(metadata$confidence_level == "strong_positive"),
      sum(metadata$confidence_level == "weak_positive")
    )
  )

  qc_final_filter_losses <- metadata_before_filter %>%
    count(class, type_group, label, confidence_level, final_filter_reason, sort = TRUE)

  qc_by_label_after <- metadata %>%
    group_by(class, type_group, label, confidence_level) %>%
    summarise(
      n_entries = n(),
      n_unique_sequences = n_distinct(sequence_md5),
      min_length = min(length, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE),
      max_length = max(length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(class, type_group, label)

  write.table(qc_selection_counts, file.path(qcdir, "qc_01_selection_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(qc_final_filter_losses, file.path(qcdir, "qc_02_final_filter_losses.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(qc_by_label_after, file.path(qcdir, "qc_03_by_label_after_filter.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\nCounts per selection step:\n")
  print(qc_selection_counts, n = Inf)
  cat("\nBy label after the final filter:\n")
  print(qc_by_label_after, n = Inf)


  # ----------------------------------------------------------
  # 5.9. Export
  # ----------------------------------------------------------

  write.table(metadata, file.path(outdir, "metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  write_fasta_final(metadata, file.path(outdir, "pos.fasta"))

  summary_host <- tibble(
    host_simple = HOST,
    n_assays_raw = nrow(df_host),
    n_df_base = nrow(df_base),
    n_positives = nrow(metadata),
    n_strong = sum(metadata$confidence_level == "strong_positive"),
    n_weak   = sum(metadata$confidence_level == "weak_positive"),
    n_epitope        = sum(metadata$type_group == "epitope"),
    n_fragment       = sum(metadata$type_group == "fragment"),
    n_source_antigen = sum(metadata$type_group == "source_antigen"),
    n_assay_antigen  = sum(metadata$type_group == "assay_antigen")
  )

  write.table(summary_host, file.path(outdir, "summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\nHost summary:\n")
  print(summary_host)
}

cat("\nOK. Positives mined into:", OUTDIR_ROOT, "/\n")
