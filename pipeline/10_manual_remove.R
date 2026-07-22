#!/usr/bin/env Rscript
# 10_manual_remove.R
# Standalone R script to remove manually flagged sequences from the final
# alignment and its taxonomy.
#
# Inputs (relative to --base-dir):
#   manual_removal.txt
#   alignment_final_noguide.fasta
#   new_taxonomy.tsv
#
# Outputs (in 10_manual_removed/):
#   alignment_final_cleaned.fasta
#   taxonomy_final_cleaned.tsv
#
# Usage:
#   Rscript 10_manual_remove.R [--base-dir <path>]
#
# Required package: phylotools

if (!requireNamespace("phylotools", quietly = TRUE)) {
  stop("The 'phylotools' R package is required.\n",
       "Install it with:  R -e \"install.packages('phylotools', repos='https://cloud.r-project.org')\"")
}
suppressPackageStartupMessages(library(phylotools))

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

removal_file <- file.path(base_dir, "manual_removal.txt")
fasta_in     <- file.path(base_dir, "alignment_final.fasta")
tax_in       <- file.path(base_dir, "new_taxonomy.tsv")

out_dir <- file.path(base_dir, "10_manual_removed")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fasta_out <- file.path(out_dir, "alignment_final_cleaned.fasta")
tax_out   <- file.path(out_dir, "taxonomy_final_cleaned.tsv")

cat(sprintf("[%s] Step 10: manual removal of flagged sequences\n", format(Sys.time())))
cat(sprintf("  Base directory: %s\n", normalizePath(base_dir)))
cat(sprintf("  Removal list : %s\n", removal_file))
cat(sprintf("  Input FASTA  : %s\n", fasta_in))
cat(sprintf("  Input tax    : %s\n", tax_in))

# ---- Read removal IDs --------------------------------------------------------
removal_lines <- readLines(removal_file)
removal_ids <- unique(trimws(removal_lines))
# Drop blank lines and comment lines
removal_ids <- removal_ids[nzchar(removal_ids) & !grepl("^#", removal_ids)]
cat(sprintf("  %d unique IDs to remove\n", length(removal_ids)))

# ---- Filter FASTA ------------------------------------------------------------
cat(sprintf("[%s] Reading FASTA ...\n", format(Sys.time())))
seqs <- read.fasta(fasta_in)
cat(sprintf("  Input sequences: %d\n", nrow(seqs)))

seqs_keep <- seqs[!(seqs$seq.name %in% removal_ids), ]
cat(sprintf("  After removal  : %d\n", nrow(seqs_keep)))

dat2fasta(seqs_keep, outfile = fasta_out)
cat(sprintf("  Wrote: %s\n", fasta_out))

# ---- Filter taxonomy ---------------------------------------------------------
cat(sprintf("[%s] Reading taxonomy ...\n", format(Sys.time())))
tax <- read.delim(tax_in, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cat(sprintf("  Input taxa: %d\n", nrow(tax)))

# First column is the Feature ID regardless of its name
feature_col <- names(tax)[1]
tax_keep <- tax[!(tax[[feature_col]] %in% removal_ids), ]
cat(sprintf("  After removal: %d\n", nrow(tax_keep)))

write.table(tax_keep, file = tax_out, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", tax_out))

# ---- Sanity check ------------------------------------------------------------
fasta_ids <- seqs_keep$seq.name
tax_ids <- tax_keep[[feature_col]]
cat(sprintf("[%s] Summary:\n", format(Sys.time())))
cat(sprintf("  Sequences in FASTA but not taxonomy: %d\n",
            length(setdiff(fasta_ids, tax_ids))))
cat(sprintf("  Taxonomy IDs but not in FASTA      : %d\n",
            length(setdiff(tax_ids, fasta_ids))))
cat(sprintf("  Common IDs                         : %d\n",
            length(intersect(fasta_ids, tax_ids))))

cat(sprintf("[%s] Step 10 complete. Results in %s\n",
            format(Sys.time()), normalizePath(out_dir)))
