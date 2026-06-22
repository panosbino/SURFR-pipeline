# =============================================================================
# SURFR Pipeline - Step 4: FindCancerSpecificRNAs.R
#
# Identifies cancer-specific k-mers by:
#   1. Loading merged k-mer count tables from TCGA and CPTAC (Step 3 outputs)
#   2. Filtering for cancer-enriched k-mers (cancer > 200, enrichment > 40,
#      adjacent < 100) independently in each cohort
#   3. Intersecting the filtered k-mers across TCGA and CPTAC
#   4. Loading SRA (non-cancer control) k-mer counts and removing any
#      intersected k-mer with SRA counts >= 200
#   5. Saving filtered results and merging overlapping k-mers with dekupl-mergeTags
#
# Usage (Rscript):
#   Rscript FindCancerSpecificRNAs.R \
#     <project_id> \
#     <merged_tables_dir> \
#     <sra_kmer_table> \
#     <metadata_dir> \
#     <output_dir> \
#     <dekupl_mergetags_path>
#
# Arguments:
#   project_id           - Project identifier, e.g. LUAD
#   merged_tables_dir    - Directory containing the Step 3 merged k-mer tables
#                          (expects files named all_<project>_TCGA_17mers_merged.txt
#                           and all_<project>_CPTAC_17mers_merged.txt)
#   sra_kmer_table       - Path to the SRA non-cancer k-mer count table
#                          (tab-delimited, columns: kmer, sra_count)
#   metadata_dir         - Directory containing TCGA and CPTAC metadata files
#   output_dir           - Directory for all output files and plots
#   dekupl_mergetags_path - Full path to the dekupl-mergeTags executable
# =============================================================================

# ----------------------------
# LOAD LIBRARIES
# ----------------------------
library(tidyverse)   # data wrangling, ggplot2, dplyr, etc.
library(paletteer)   # additional color palettes
library(arrow)       # write_feather output format
library(ggvenn)      # Venn diagrams

# ----------------------------
# PARSE COMMAND-LINE ARGUMENTS
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "ERROR: Expected 6 arguments.",
    "Usage: Rscript FindCancerSpecificRNAs.R",
    "  <project_id> <merged_tables_dir> <sra_kmer_table>",
    "  <metadata_dir> <output_dir> <dekupl_mergetags_path>"
  ))
}

project              <- args[1]
merged_tables_dir    <- args[2]
sra_kmer_table_path  <- args[3]
metadata_dir         <- args[4]
output_dir           <- args[5]
dekupl_path          <- args[6]

kmer_length <- 17   # FIX: corrected typo from 'kmer_lenght'

# ----------------------------
# CREATE OUTPUT DIRECTORIES
# ----------------------------
analysis_dir <- file.path(output_dir, "Analysis")
plots_dir    <- file.path(output_dir, "Plots")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,    recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[%s] Output directory: %s\n", Sys.time(), output_dir))

# ----------------------------
# HELPER: POINT DENSITY
# (replaces get_density() from Utils_new.R for portability)
# ----------------------------
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x, y, n = n,
                      lims = c(range(x), range(y)))
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Custom color palette for density plots
colors <- c("#3d3d3d","#747474","#989898","#cacaca","#c8c8c8","#dadada","#e4e4e4")

# ----------------------------
# LOAD METADATA
# ----------------------------
cat(sprintf("[%s] Loading metadata...\n", Sys.time()))

metadata_TCGA <- read.delim(
  file.path(metadata_dir, sprintf("clean_metadata_long_TCGA_%s.txt", project)),
  sep = ";"
)

metadata_CPTAC <- read.delim(
  file.path(metadata_dir, sprintf("clean_metadata_long_CPTAC_%s.txt", project)),
  sep = " "
)

# Compute sample ratios (cancer / adjacent-normal) for reference
sample_ratio_CPTAC <- metadata_CPTAC |> pull(Sample_Type) |> table()
sample_ratio_CPTAC <- sample_ratio_CPTAC[1] / sample_ratio_CPTAC[2]

sample_ratio_TCGA  <- metadata_TCGA  |> pull(Sample_Type) |> table()
sample_ratio_TCGA  <- sample_ratio_TCGA[1] / sample_ratio_TCGA[2]

# ----------------------------
# LOAD K-MER DATA — TCGA
# ----------------------------
cat(sprintf("[%s] Loading TCGA k-mer table...\n", Sys.time()))

tcga_file <- file.path(merged_tables_dir,
                       sprintf("all_%s_TCGA_%dmers_merged.txt", project, kmer_length))
if (!file.exists(tcga_file)) stop(paste("TCGA k-mer file not found:", tcga_file))

df_TCGA <- read.delim(tcga_file, header = FALSE)
colnames(df_TCGA) <- c("kmers", "sums_TCGA_Cancer_KMC", "sums_TCGA_Healthy_KMC")
df_TCGA[is.na(df_TCGA)] <- 0

df_TCGA <- df_TCGA |>
  mutate(enrichment_TCGA = sums_TCGA_Cancer_KMC / sums_TCGA_Healthy_KMC)

# Apply enrichment filters
tcga_filter <- df_TCGA |>
  filter(sums_TCGA_Cancer_KMC > 200) |>
  filter(enrichment_TCGA > 40) |>
  filter(sums_TCGA_Healthy_KMC < 100)

cat(sprintf("[%s] TCGA: %d k-mers passed enrichment filters\n",
            Sys.time(), nrow(tcga_filter)))

# ----------------------------
# LOAD K-MER DATA — CPTAC
# ----------------------------
cat(sprintf("[%s] Loading CPTAC k-mer table...\n", Sys.time()))

cptac_file <- file.path(merged_tables_dir,
                        sprintf("all_%s_CPTAC_%dmers_merged.txt", project, kmer_length))
if (!file.exists(cptac_file)) stop(paste("CPTAC k-mer file not found:", cptac_file))

df_CPTAC <- read.delim(cptac_file, header = FALSE)
colnames(df_CPTAC) <- c("kmers", "sums_CPTAC_Cancer_KMC", "sums_CPTAC_Healthy_KMC")
df_CPTAC[is.na(df_CPTAC)] <- 0

df_CPTAC <- df_CPTAC |>
  mutate(enrichment_CPTAC = sums_CPTAC_Cancer_KMC / sums_CPTAC_Healthy_KMC)

# Apply enrichment filters
cptac_filter <- df_CPTAC |>
  filter(sums_CPTAC_Cancer_KMC > 200) |>
  filter(enrichment_CPTAC > 40) |>
  filter(sums_CPTAC_Healthy_KMC < 100)

cat(sprintf("[%s] CPTAC: %d k-mers passed enrichment filters\n",
            Sys.time(), nrow(cptac_filter)))

# ----------------------------
# INTERSECT FILTERED K-MERS ACROSS COHORTS
# ----------------------------
cat(sprintf("[%s] Intersecting TCGA and CPTAC filtered k-mers...\n", Sys.time()))

intersected_kmers <- intersect(cptac_filter$kmers, tcga_filter$kmers)
cat(sprintf("[%s] %d k-mers shared between TCGA and CPTAC cohorts\n",
            Sys.time(), length(intersected_kmers)))

# FIX: Build a merged df joining TCGA and CPTAC counts — 'df' was previously
# undefined, causing a runtime crash. Full join ensures all intersected k-mers
# are present regardless of which table they originate from.
df <- dplyr::full_join(df_TCGA, df_CPTAC, by = "kmers")
df[is.na(df)] <- 0

intersected_df <- df |> dplyr::filter(kmers %in% intersected_kmers)

# ----------------------------
# LOAD SRA (NON-CANCER CONTROL) K-MER COUNTS
# FIX: SRA data was referenced but never loaded in original script
# ----------------------------
cat(sprintf("[%s] Loading SRA k-mer counts...\n", Sys.time()))

if (!file.exists(sra_kmer_table_path)) {
  stop(paste("SRA k-mer table not found:", sra_kmer_table_path))
}

sra_df <- read.delim(sra_kmer_table_path, header = FALSE)
colnames(sra_df) <- c("kmers", "sums_SRA")

# Join SRA counts onto intersected k-mers; absent k-mers get count 0
intersected_df <- intersected_df |>
  dplyr::left_join(sra_df, by = "kmers") |>
  tidyr::replace_na(list(sums_SRA = 0))

# ----------------------------
# FILTER USING SRA COUNTS
# Remove k-mers that are abundant in non-cancer controls
# ----------------------------
sra_cutoff <- 200

cancer_specific_df <- intersected_df |>
  filter(sums_SRA < sra_cutoff)

cat(sprintf("[%s] %d k-mers remain after SRA filtering (cutoff: %d counts)\n",
            Sys.time(), nrow(cancer_specific_df), sra_cutoff))

# Save intersected k-mer list
write(
  x    = intersected_kmers,
  file = file.path(analysis_dir, sprintf("intersected_kmers_%s.txt", project))
)

# ----------------------------
# PLOTS: SAMPLE DISTRIBUTIONS
# ----------------------------
cat(sprintf("[%s] Generating sample distribution plots...\n", Sys.time()))

bar_TCGA_H_C <- ggplot() +
  geom_bar(data    = metadata_TCGA,
           mapping = aes(y = Sample_Type, fill = Sample_Type),
           width = 0.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#ce3701", "#72cacf")) +
  theme(
    axis.title.x       = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y       = element_text(size = 20, margin = margin(r = 20)),
    axis.text          = element_text(size = 18, colour = "black"),
    axis.ticks.y       = element_line(linewidth = 1),
    axis.ticks.x       = element_line(linewidth = 1),
    plot.title         = element_text(size = 18),
    axis.ticks.length  = unit(.2, "cm"),
    axis.title.y.left  = element_blank(),
    axis.title.x.bottom = element_blank(),
    legend.text        = element_text(size = 18),
    legend.title       = element_blank()
  ) +
  ggtitle("Number of TCGA samples")

ggsave(plot     = bar_TCGA_H_C,
       path     = plots_dir,
       filename = sprintf("barplot_%s_TCGA_Healthy_Cancer.png", project),
       dpi = 600, device = "png", width = 8, height = 4)
ggsave(plot     = bar_TCGA_H_C,
       path     = plots_dir,
       filename = sprintf("barplot_%s_TCGA_Healthy_Cancer.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 2)

bar_CPTAC_H_C <- ggplot() +
  geom_bar(data    = metadata_CPTAC,
           mapping = aes(y = Sample_Type, fill = Sample_Type),
           width = 0.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#ce3701", "#72cacf")) +
  theme(
    axis.title.x        = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y        = element_text(size = 20, margin = margin(r = 20)),
    axis.text           = element_text(size = 18, colour = "black"),
    axis.ticks.y        = element_line(linewidth = 1),
    axis.ticks.x        = element_line(linewidth = 1),
    plot.title          = element_text(size = 18),
    axis.ticks.length   = unit(.2, "cm"),
    axis.title.y.left   = element_blank(),
    axis.title.x.bottom = element_blank(),
    legend.text         = element_text(size = 18),
    legend.title        = element_blank()
  ) +
  ggtitle("Number of CPTAC samples")

ggsave(plot     = bar_CPTAC_H_C,
       path     = plots_dir,
       filename = sprintf("barplot_%s_CPTAC_Healthy_Cancer.png", project),
       dpi = 600, device = "png", width = 8, height = 4)
ggsave(plot     = bar_CPTAC_H_C,
       path     = plots_dir,
       filename = sprintf("barplot_%s_CPTAC_Healthy_Cancer.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 2)

# ----------------------------
# PLOTS: SCATTERPLOTS — TCGA
# ----------------------------
cat(sprintf("[%s] Generating TCGA scatterplot...\n", Sys.time()))

df_TCGA$tcga_filter <- ifelse(df_TCGA$kmers %in% tcga_filter$kmers, "passed", "failed")
df_TCGA$density_TCGA <- get_density(
  log10(df_TCGA$sums_TCGA_Cancer_KMC + 1),
  log10(df_TCGA$sums_TCGA_Healthy_KMC + 1)
)

p_all_tcga_filt <- ggplot() +
  geom_point(data    = df_TCGA,
             mapping = aes(x     = log10(sums_TCGA_Healthy_KMC + 1),
                           y     = log10(sums_TCGA_Cancer_KMC + 1),
                           fill  = log10(density_TCGA)),
             size = 2, shape = 21, stroke = NA) +
  geom_point(data    = tcga_filter,
             mapping = aes(x = log10(sums_TCGA_Healthy_KMC + 1),
                           y = log10(sums_TCGA_Cancer_KMC + 1)),
             size = 2, shape = 21, fill = "#febf38", stroke = NA) +
  scale_x_continuous(limits = c(0, log10(max(df_TCGA$sums_TCGA_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  scale_y_continuous(limits = c(0, log10(max(df_TCGA$sums_TCGA_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  theme(
    panel.border      = element_blank(),
    panel.grid        = element_blank(),
    axis.title        = element_blank(),
    axis.text         = element_text(size = 18, colour = "black"),
    axis.ticks.y      = element_line(linewidth = 1),
    axis.ticks.x      = element_line(linewidth = 1),
    plot.title        = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    panel.background  = element_rect(fill = "transparent"),
    plot.background   = element_rect(fill = "transparent", color = NA),
    legend.position   = "none"
  ) +
  scale_fill_gradientn(colours = rev(colors))

ggsave(plot     = p_all_tcga_filt,
       path     = plots_dir,
       filename = sprintf("scatterplot_%s_TCGA_FILT_all.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 6, bg = "transparent")

# ----------------------------
# PLOTS: SCATTERPLOTS — CPTAC
# ----------------------------
cat(sprintf("[%s] Generating CPTAC scatterplot...\n", Sys.time()))

df_CPTAC$cptac_filter <- ifelse(df_CPTAC$kmers %in% cptac_filter$kmers, "passed", "failed")
df_CPTAC$density_CPTAC <- get_density(
  log10(df_CPTAC$sums_CPTAC_Cancer_KMC + 1),
  log10(df_CPTAC$sums_CPTAC_Healthy_KMC + 1)
)

p_all_CPTAC_filt <- ggplot() +
  geom_point(data    = df_CPTAC,
             mapping = aes(x    = log10(sums_CPTAC_Healthy_KMC + 1),
                           y    = log10(sums_CPTAC_Cancer_KMC + 1),
                           fill = log10(density_CPTAC)),
             size = 2, shape = 21, stroke = NA) +
  geom_point(data    = cptac_filter,
             mapping = aes(x = log10(sums_CPTAC_Healthy_KMC + 1),
                           y = log10(sums_CPTAC_Cancer_KMC + 1)),
             size = 2, shape = 21, fill = "#febf38", stroke = NA) +
  scale_x_continuous(limits = c(0, log10(max(df_CPTAC$sums_CPTAC_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  scale_y_continuous(limits = c(0, log10(max(df_CPTAC$sums_CPTAC_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  theme(
    panel.border      = element_blank(),
    panel.grid        = element_blank(),
    axis.title        = element_blank(),
    axis.text         = element_text(size = 18, colour = "black"),
    axis.ticks.y      = element_line(linewidth = 1),
    axis.ticks.x      = element_line(linewidth = 1),
    plot.title        = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    panel.background  = element_rect(fill = "transparent"),
    plot.background   = element_rect(fill = "transparent", color = NA),
    legend.position   = "none"
  ) +
  scale_fill_gradientn(colours = rev(colors))

ggsave(plot     = p_all_CPTAC_filt,
       path     = plots_dir,
       filename = sprintf("scatterplot_%s_CPTAC_FILT_all.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 6, bg = "transparent")

# ----------------------------
# PLOTS: VENN DIAGRAM
# ----------------------------
cat(sprintf("[%s] Generating Venn diagram...\n", Sys.time()))

venn_list <- list(
  CPTAC = cptac_filter$kmers,
  TCGA  = tcga_filter$kmers
)

ven <- ggvenn(data = venn_list, fill_color = c("#40af75", "#1e6bb0"))

ggsave(plot     = ven,
       path     = plots_dir,
       filename = sprintf("venn_%s_CPTAC_TCGA.png", project),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot     = ven,
       path     = plots_dir,
       filename = sprintf("venn_%s_CPTAC_TCGA.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# PLOTS: SRA DENSITY PLOT
# ----------------------------
cat(sprintf("[%s] Generating SRA density plot...\n", Sys.time()))

dens_plot_SRA <- ggplot() +
  geom_density(data    = intersected_df,
               mapping = aes(x = log10(sums_SRA + 1)),
               fill = "#dacdd9") +
  geom_vline(xintercept = log10(sra_cutoff),
             linetype = "dashed", linewidth = 1) +
  theme_classic() +
  theme(
    axis.title.x      = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 20)),
    axis.text         = element_text(size = 18, colour = "black"),
    axis.ticks.y      = element_line(linewidth = 1),
    axis.ticks.x      = element_line(linewidth = 1),
    plot.title        = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm")
  ) +
  xlab(label = "log10 SRA counts") +
  xlim(c(0, 1.5 * max(log10(intersected_df$sums_SRA + 1)))) +
  ggtitle(label = sprintf("%d k-mers have < %d SRA counts",
                          nrow(cancer_specific_df), sra_cutoff))

ggsave(plot     = dens_plot_SRA,
       path     = plots_dir,
       filename = sprintf("density_%s_SRA_counts.png", project),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot     = dens_plot_SRA,
       path     = plots_dir,
       filename = sprintf("density_%s_SRA_counts.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# PLOTS: VIOLIN PLOT
# ----------------------------
cat(sprintf("[%s] Generating violin plot...\n", Sys.time()))

tmp <- cancer_specific_df |>
  pivot_longer(cols = c(sums_CPTAC_Cancer_KMC,
                        sums_CPTAC_Healthy_KMC,
                        sums_TCGA_Cancer_KMC,
                        sums_TCGA_Healthy_KMC,
                        sums_SRA))

a <- tmp$name |>
  str_split_i(pattern = "sums_", i = 2) |>
  str_replace(pattern = "_KMC", replacement = "") |>
  str_replace(pattern = "_", replacement = " ") |>
  str_replace(pattern = "Healthy", replacement = "Adjacent Normal")

tmp$name <- ifelse(a == "SRA", "Non-Cancer SRA", a)
tmp$name <- factor(tmp$name,
                   levels = c("CPTAC Cancer", "TCGA Cancer",
                              "CPTAC Adjacent Normal", "TCGA Adjacent Normal",
                              "Non-Cancer SRA"))

tmp$condition <- ifelse(tmp$name %in% c("CPTAC Cancer", "TCGA Cancer"),
                        "cancer", "non-cancer")

violin_counts <- ggplot() +
  geom_violin(data    = tmp,
              mapping = aes(x = name, y = log10(value + 1), fill = condition)) +
  geom_jitter(data    = tmp,
              mapping = aes(x = name, y = log10(value + 1)),
              size = 0.2, position = position_jitter(0.2)) +
  theme_bw() +
  theme(
    panel.background    = element_blank(),
    panel.border        = element_rect(fill = NA, linewidth = 2),
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    axis.title.x        = element_blank(),
    axis.title.y        = element_text(size = 20, margin = margin(r = 20)),
    axis.text           = element_text(size = 18, colour = "black"),
    axis.ticks.y        = element_line(linewidth = 1),
    axis.ticks.x        = element_line(linewidth = 1),
    plot.title          = element_text(size = 18),
    axis.ticks.length   = unit(.2, "cm"),
    legend.title        = element_blank(),
    legend.text         = element_text(size = 18),
    axis.text.x         = element_text(angle = 45, vjust = 0.7)
  ) +
  scale_fill_manual(values = c("#ce3701", "#72cacf")) +
  ylab(label = "log10 counts") +
  scale_x_discrete(labels = function(x) { sub("\\s", "\n", x) })

ggsave(plot     = violin_counts,
       path     = plots_dir,
       filename = sprintf("violin_%s_candidate_counts.png", project),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot     = violin_counts,
       path     = plots_dir,
       filename = sprintf("violin_%s_candidate_counts.pdf", project),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# SAVE FILTERED RESULTS
# ----------------------------
cat(sprintf("[%s] Saving filtered results...\n", Sys.time()))

arrow::write_feather(
  cancer_specific_df,
  sink = file.path(analysis_dir,
                   sprintf("cancer_enriched_filtered_%dmers_%s.arrow",
                           kmer_length, project))
)

write_delim(
  x    = cancer_specific_df,
  file = file.path(analysis_dir,
                   sprintf("cancer_enriched_filtered_%dmers_%s.tsv",
                           kmer_length, project)),
  delim = "\t"
)

# ----------------------------
# MERGE K-MERS WITH dekupl-mergeTags
# FIX: path is now passed as a command-line argument, not hardcoded
# ----------------------------
cat(sprintf("[%s] Running dekupl-mergeTags...\n", Sys.time()))

if (!file.exists(dekupl_path)) {
  stop(paste("dekupl-mergeTags executable not found:", dekupl_path))
}

input_tsv  <- file.path(analysis_dir,
                        sprintf("cancer_enriched_filtered_%dmers_%s.tsv",
                                kmer_length, project))
output_tsv <- file.path(analysis_dir,
                        sprintf("cancer_enriched_filtered_sequences_%s.tsv", project))

system(sprintf("%s -k %d -m 8 -n %s > %s",
               dekupl_path, kmer_length, input_tsv, output_tsv))

# Load merged candidate sequences
cancer_specific_sequences <- read.delim(output_tsv)

cat(sprintf("[%s] %d candidate sequences identified after merging.\n",
            Sys.time(), nrow(cancer_specific_sequences)))

cat(sprintf("[%s] Step 4 (FindCancerSpecificRNAs) complete.\n", Sys.time()))
cat(sprintf("[%s] Final output: %s\n", Sys.time(), output_tsv))
