# =============================================================================
# SURFR v1 — regression test against the published LUAD result
#
# Runs the non-cancer (SRA) filter and the dekupl-mergeTags merge exactly as
# FindCancerSpecificRNAs.r does (the script's own expressions are evaluated,
# nothing is re-implemented), starting from the 156 LUAD k-mers that passed the
# SRA filter in the published analysis. It checks that the output is exactly
# the 73 published LUAD sequences (Supplementary Table 2; one of them is listed
# there under LUSC as oncRNA-8) with the same representative k-mers.
#
# The cohort count columns are placeholders: they do not influence which
# sequences are produced or which k-mer represents each sequence.
#
# Usage (from the repository root, inside the container):
#   Rscript tests/test_step4_LUAD_published.R /opt/dekupl/bin/mergeTags
# Exit status 0 = pass, 1 = fail.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args        <- commandArgs(trailingOnly = TRUE)
dekupl_path <- if (length(args) >= 1) args[1] else "/opt/dekupl/bin/mergeTags"
script      <- "FindCancerSpecificRNAs.r"
fixtures    <- "tests/data"
if (!file.exists(script)) stop("Run this test from the repository root.")

fail <- function(...) { cat("FAIL:", ..., "\n"); quit(status = 1) }

# --- Inputs ------------------------------------------------------------------
published <- read.delim(file.path(fixtures, "published_LUAD_sequences.tsv"),
                        colClasses = "character")
sra_pub   <- read.table(file.path(fixtures, "sums_SRA_filtered_LUAD.txt"),
                        col.names = c("kmers", "sums_SRA"),
                        colClasses = c("character", "numeric"))

# Two extra k-mers that must be REMOVED by the SRA filter (boundary: >= 200)
extra <- data.frame(kmers    = c("ACGTACGTACGTACGTA", "TTTTTCCCCCGGGGGAA"),
                    sums_SRA = c(200, 5000))

work <- tempfile("surfr_test_")
dir.create(work)
sra_kmer_table_path <- file.path(work, "sra_counts.txt")
# space-separated, like the original table
write.table(rbind(sra_pub, extra), sra_kmer_table_path, sep = " ",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

all_kmers <- sort(c(sra_pub$kmers, extra$kmers))
intersected_df <- data.frame(
  kmers                  = all_kmers,
  sums_TCGA_Cancer_KMC   = 1000, sums_TCGA_Healthy_KMC  = 5,
  enrichment_TCGA        = 200,
  sums_CPTAC_Cancer_KMC  = 800,  sums_CPTAC_Healthy_KMC = 0,
  enrichment_CPTAC       = Inf
)
analysis_dir <- work
project      <- "LUAD"
kmer_length  <- 17

# --- Evaluate the script's own SRA-filter and merge expressions --------------
ex  <- parse(script, keep.source = TRUE)
src <- vapply(attr(ex, "srcref"),
              function(s) paste(as.character(s), collapse = "\n"), "")
pattern <- paste0("^(sra_df <-|intersected_df <- intersected_df|sra_cutoff <-|",
                  "cancer_specific_df <-|write_delim|input_tsv +<-|output_tsv <-|",
                  "mergetags_|if \\(mergetags_status|cancer_specific_sequences <-)")
use <- grepl(pattern, src)
if (sum(use) != 12) {
  fail(sprintf(paste("expected to find 12 SRA-filter/merge expressions in %s,",
                     "found %d. The script structure changed; update this test."),
               script, sum(use)))
}
for (i in which(use)) eval(ex[[i]])

# --- Checks ------------------------------------------------------------------
if (nrow(cancer_specific_df) != nrow(sra_pub))
  fail(sprintf("%d k-mers passed the SRA filter, expected %d",
               nrow(cancer_specific_df), nrow(sra_pub)))
if (any(extra$kmers %in% cancer_specific_df$kmers))
  fail("k-mers with SRA counts >= 200 were not removed")

got <- setNames(cancer_specific_sequences$kmers, cancer_specific_sequences$contig)
exp <- setNames(published$representative_17mer, published$contig)

missing    <- setdiff(names(exp), names(got))
unexpected <- setdiff(names(got), names(exp))
if (length(missing) || length(unexpected))
  fail(sprintf("sequence sets differ. Missing: %s | Unexpected: %s",
               paste(missing, collapse = ","), paste(unexpected, collapse = ",")))

wrong_rep <- names(exp)[got[names(exp)] != exp]
if (length(wrong_rep))
  fail(sprintf("representative k-mer differs for: %s", paste(wrong_rep, collapse = ",")))

cat(sprintf("PASS: %d sequences and representative k-mers identical to the published LUAD result\n",
            length(exp)))
