#!/usr/bin/env Rscript
# 13_taxonomy_test.R
# Standalone R script to test the short- and long-amplicon dereplicated
# databases as DADA2 taxonomy references.
#
# For each region:
#   1. Build a DADA2-formatted reference FASTA from dereplicated sequences
#      and their taxonomy strings.
#   2. Classify all query sequences with DADA2::assignTaxonomy.
#   3. Merge classification results with original query IDs and write CSV.
#
# Usage:
#   Rscript 13_taxonomy_test.R [--base-dir <path>]
#
# Required packages: phylotools, dada2, dplyr

suppressPackageStartupMessages({
  library(phylotools)
  library(dada2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

# Default base directory is the current working directory
base_dir <- "."

parse_arg <- function(args, flag, default) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

if ("--base-dir" %in% args) {
  base_dir <- parse_arg(args, "--base-dir", base_dir)
}

out_dir <- file.path(base_dir, "13_taxonomy_test")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat(sprintf("[%s] Step 13: testing amplicon databases as DADA2 taxonomy references\n",
            format(Sys.time())))
cat(sprintf("  Base directory: %s\n", normalizePath(base_dir)))

# ---- Helper function ---------------------------------------------------------

test_region <- function(region) {
  cat(sprintf("\n[%s] --- %s amplicon ---\n", format(Sys.time()), region))

  region_dir <- file.path(base_dir, "12_amplicon_derep", region)
  queries   <- file.path(region_dir, sprintf("seqs_%s.fasta", region))
  ref_seq   <- file.path(region_dir, sprintf("seqs_%s_derep.fasta", region))
  # Ultimate expected taxonomy (full dereplicated database)
  ref_tax   <- file.path(base_dir, "taxonomy_derep_final.tsv")

  cat(sprintf("  Queries:  %s\n", queries))
  cat(sprintf("  Ref seqs: %s\n", ref_seq))
  cat(sprintf("  Ref tax:  %s\n", ref_tax))

  # Read the ultimate reference taxonomy (one row per original reference ID)
  tax_df <- read.delim(ref_tax, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  colnames(tax_df) <- c("seq.name", "Taxon")
  tax_df$expected_species <- sapply(strsplit(tax_df$Taxon, ";"), tail, n = 1)

  # Read dereplicated amplicon reference sequences
  ref_seqs <- read.fasta(ref_seq)
  cat(sprintf("  Reference sequences: %d\n", nrow(ref_seqs)))

  ref_merged <- merge(ref_seqs, tax_df, by = "seq.name", all.x = FALSE)
  if (nrow(ref_merged) != nrow(ref_seqs)) {
    stop(sprintf("[%s] Some reference sequences are missing from the taxonomy file.", region))
  }

  # Write DADA2 reference FASTA (headers are taxonomy strings)
  ref_fasta <- file.path(out_dir, sprintf("seqs_%s_dada.fas", region))
  writeLines(
    paste0(">", ref_merged$Taxon, "\n", toupper(ref_merged$seq.text)),
    con = ref_fasta
  )
  cat(sprintf("  Wrote DADA2 reference: %s\n", ref_fasta))

  # Read query sequences from the trimmed (non-dereplicated) dataset
  query_seqs <- read.fasta(queries)
  cat(sprintf("  Query sequences: %d\n", nrow(query_seqs)))

  # Build a DADA2-style seqtable: one sample (row) x unique sequences (columns)
  queries_upper <- toupper(query_seqs$seq.text)
  unique_seqs <- unique(queries_upper)
  cat(sprintf("  Unique query sequences: %d\n", length(unique_seqs)))

  seq_counts <- table(factor(queries_upper, levels = unique_seqs))
  seqtable <- matrix(as.integer(seq_counts), nrow = 1)
  colnames(seqtable) <- unique_seqs
  rownames(seqtable) <- "sample1"

  # Run DADA2 assignTaxonomy
  cat(sprintf("[%s] Running DADA2 assignTaxonomy for %s ...\n",
              format(Sys.time()), region))

  tax_levels <- c("Domain", "Kingdom", "Phylum", "Class",
                  "Order", "Family", "Genus", "Species")

  tax_assigned <- assignTaxonomy(
    seqtable,
    refFasta = ref_fasta,
    taxLevels = tax_levels,
    minBoot = 0,
    tryRC = TRUE,
    outputBootstraps = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )

  # Taxonomy assignments
  tax_df_out <- as.data.frame(tax_assigned$tax)
  tax_df_out$sequence <- rownames(tax_df_out)
  rownames(tax_df_out) <- NULL

  # Bootstrap values per taxonomic level
  boot_df_out <- as.data.frame(tax_assigned$boot)
  boot_df_out$sequence <- rownames(boot_df_out)
  rownames(boot_df_out) <- NULL
  names(boot_df_out) <- c(paste0(tax_levels, "_boot"), "sequence")

  # Build per-query result. Keep one row per query sequence, with the count
  # of how many times the exact sequence appears across the query set.
  query_lookup <- query_seqs %>%
    mutate(sequence = toupper(seq.text)) %>%
    select(seq.name, sequence)

  seq_counts <- query_lookup %>%
    group_by(sequence) %>%
    summarise(count = n(), .groups = "drop")

  query_lookup <- query_lookup %>%
    left_join(seq_counts, by = "sequence") %>%
    mutate(ref_id = seq.name)

  # Look up expected species from the ultimate reference taxonomy by ID.
  # Every original reference ID has a taxonomy entry in the ultimate file.
  ref_tax_lookup <- tax_df %>%
    select(seq.name, expected_species)

  n_missing <- sum(!(query_lookup$seq.name %in% ref_tax_lookup$seq.name))
  if (n_missing > 0) {
    cat(sprintf("  Warning: %d query IDs are missing from the reference taxonomy\n", n_missing))
  }

  # Build result table
  result <- query_lookup %>%
    left_join(ref_tax_lookup, by = c("seq.name" = "seq.name")) %>%
    left_join(tax_df_out, by = "sequence") %>%
    left_join(boot_df_out, by = "sequence")

  # Species-level comparison
  result$assigned_species <- result$Species
  result$species_match <- ifelse(
    is.na(result$expected_species) | is.na(result$assigned_species),
    NA_character_,
    ifelse(result$expected_species == result$assigned_species, "yes", "no")
  )

  # Sort by sequence so identical sequences are adjacent
  result <- result %>%
    select(seq.name, sequence, count, ref_id, expected_species,
           assigned_species, species_match, everything()) %>%
    arrange(sequence)

  n_match <- sum(result$species_match == "yes", na.rm = TRUE)
  n_compare <- sum(!is.na(result$species_match))
  cat(sprintf("  Query sequences: %d\n", nrow(result)))
  cat(sprintf("  Species assignments matching expected: %d / %d\n",
              n_match, n_compare))

  # Write classification results
  out_csv <- file.path(out_dir, sprintf("test_tax_%s.csv", region))
  write.csv(result, file = out_csv, row.names = FALSE)
  cat(sprintf("  Wrote results: %s\n", out_csv))
}

# ---- Run both regions --------------------------------------------------------

test_region("short")
test_region("long")

cat(sprintf("\n[%s] Step 13 complete. Results in %s\n",
            format(Sys.time()), normalizePath(out_dir)))
cat("  test_tax_short.csv\n")
cat("  test_tax_long.csv\n")
cat("  seqs_short_dada.fas\n")
cat("  seqs_long_dada.fas\n")
