# ----------------------------
# LOAD LIBRARIES
# ----------------------------
library(tidyverse)   # data wrangling, ggplot2, dplyr, etc.
library(paletteer)   # additional color palettes for plots
library(arrow)       # efficient data format support (not directly used here, but useful)
library(ggvenn)      # Venn diagrams for overlapping sets

# ----------------------------
# SETUP
# ----------------------------
setwd("~/Desktop/")  # set working directory

# Custom color palette (unused in this script, but available)
colors <- c("#3d3d3d","#747474","#989898","#cacaca","#c8c8c8","#dadada","#e4e4e4")

# Load utility functions from external script
source("cancer_kmers/scripts/Utils_new.R")

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
args <- commandArgs(trailingOnly=TRUE)  # parse command line arguments
project <- args[1]                      # project name (e.g. "LUAD")
kmer_lenght <- 17                       # fixed k-mer length for analysis

# Define directories
main.dir <- paste0("~/Desktop/cancer_kmers/",project,"_kmers/")
analysis.dir <- paste0(main.dir, "/Analysis")
plots.dir <- "~/Desktop/Cancer_paper/Plots_new_method/"

# Create directories if they don't exist
system(paste0('mkdir -p ', analysis.dir))
system(paste0('mkdir -p ', plots.dir))

# ----------------------------
# LOAD METADATA
# ----------------------------
# TCGA and CPTAC metadata contain information about sample types (cancer/healthy)
metadata_TCGA <- read.delim(
  paste0(main.dir,"TCGA/metadata_TCGA_",project,"/clean_metadata_long_TCGA_",project,".txt"), 
  sep = ";"
)
metadata_CPTAC <- read.delim(
  paste0(main.dir,"CPTAC/metadata_CPTAC_",project,"/clean_metadata_long_CPTAC_",project,".txt"), 
  sep = " "
)

# Compute ratios of cancer vs healthy samples in both datasets
sample_ratio_CPTAC <- metadata_CPTAC |> pull(Sample_Type) |> table()
sample_ratio_CPTAC <- sample_ratio_CPTAC[1]/sample_ratio_CPTAC[2]

sample_ratio_TCGA <- metadata_TCGA |> pull(Sample_Type) |> table()
sample_ratio_TCGA <- sample_ratio_TCGA[1]/sample_ratio_TCGA[2]

# ----------------------------
# LOAD K-MER DATA (TCGA)
# ----------------------------
df_TCGA <- read.delim(
  paste0("~/Desktop/all_merged_tables/all_",project,"_TCGA_17mers_merged.txt"), 
  header = FALSE
)
colnames(df_TCGA) <- c("kmers","sums_TCGA_Cancer_KMC","sums_TCGA_Healthy_KMC")
df_TCGA[is.na(df_TCGA)] <- 0   # replace NAs with 0

# Compute enrichment score (cancer count / healthy count)
df_TCGA <- df_TCGA |>
  mutate(enrichment_TCGA = sums_TCGA_Cancer_KMC / sums_TCGA_Healthy_KMC)

# Apply filters to select enriched k-mers
tcga_filter <- df_TCGA |> 
  filter(sums_TCGA_Cancer_KMC > 200) |>    # cancer signal must be at least 200 counts
  filter(enrichment_TCGA > 40) |>          # enrichment must be high
  filter(sums_TCGA_Healthy_KMC < 100)      # avoid high counts in healthy samples

# ----------------------------
# PLOT SAMPLE DISTRIBUTIONS (TCGA)
# ----------------------------
bar_TCGA_H_C <- ggplot() +
  geom_bar(data = metadata_TCGA, mapping = aes(y = Sample_Type, fill = Sample_Type),
           width = 0.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#ce3701","#72cacf")) +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20)),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    axis.title.y.left = element_blank(),
    axis.title.x.bottom = element_blank(),
    legend.text = element_text(size =18),
    legend.title = element_blank()
  ) +
  ggtitle("Number of TCGA samples")

# Save TCGA sample distribution plots
ggsave(plot = bar_TCGA_H_C, path = plots.dir, 
       filename = paste0("barplot_", project, "_TCGA_Healthy_Cancer.png"),
       dpi = 600, device = "png", width = 8, height = 4)
ggsave(plot = bar_TCGA_H_C, path = plots.dir, 
       filename = paste0("barplot_", project, "_TCGA_Healthy_Cancer.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 2)

# ----------------------------
# PLOT SAMPLE DISTRIBUTIONS (CPTAC)
# ----------------------------
bar_CPTAC_H_C <- ggplot() +
  geom_bar(data = metadata_CPTAC, mapping = aes(y = Sample_Type, fill = Sample_Type),
           width = 0.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#ce3701","#72cacf")) +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20)),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    axis.title.y.left = element_blank(),
    axis.title.x.bottom = element_blank(),
    legend.text = element_text(size =18),
    legend.title = element_blank()
  ) +
  ggtitle("Number of CPTAC samples")

# Save CPTAC sample distribution plots
ggsave(plot = bar_CPTAC_H_C, path = plots.dir, 
       filename = paste0("barplot_", project, "_CPTAC_Healthy_Cancer.png"),
       dpi = 600, device = "png", width = 8, height = 4)
ggsave(plot = bar_CPTAC_H_C, path = plots.dir, 
       filename = paste0("barplot_", project, "_CPTAC_Healthy_Cancer.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 2)

# ----------------------------
# LOAD K-MER DATA (CPTAC)
# ----------------------------
df_CPTAC <- read.delim(
  paste0("~/Desktop/all_merged_tables/all_",project,"_CPTAC_17mers_merged.txt"), 
  header = FALSE
)
colnames(df_CPTAC) <- c("kmers","sums_CPTAC_Cancer_KMC","sums_CPTAC_Healthy_KMC")
df_CPTAC[is.na(df_CPTAC)] <- 0

# Compute enrichment score
df_CPTAC <- df_CPTAC |>
  mutate(enrichment_CPTAC = sums_CPTAC_Cancer_KMC / sums_CPTAC_Healthy_KMC)

# Apply filters
cptac_filter <- df_CPTAC |> 
  filter(sums_CPTAC_Cancer_KMC > 200) |>
  filter(enrichment_CPTAC > 40) |>
  filter(sums_CPTAC_Healthy_KMC < 100)


# ----------------------------
# SCATTERPLOTS FOR CANCER VS ADJACENT
# ----------------------------

# ---- TCGA ----
# Mark k-mers that pass TCGA filtering
df_TCGA$tcga_filter <- ifelse(df_TCGA$kmers %in% tcga_filter$kmers, "passed", "failed")

# Calculate point density for scatterplot visualization
df_TCGA$density_TCGA <- get_density(
  log10(df_TCGA$sums_TCGA_Cancer_KMC + 1),
  log10(df_TCGA$sums_TCGA_Healthy_KMC + 1)
)

# Scatterplot: all TCGA vs filtered highlighted
p_all_tcga_filt <- ggplot() +
  geom_point_rast(data = df_TCGA,
                  mapping = aes(x = log10(sums_TCGA_Healthy_KMC + 1),
                                y = log10(sums_TCGA_Cancer_KMC + 1),
                                fill = log10(density_TCGA)),
                  size = 2, shape = 21, stroke = NA) +
  geom_point_rast(data = tcga_filter,
                  mapping = aes(x = log10(sums_TCGA_Healthy_KMC + 1),
                                y = log10(sums_TCGA_Cancer_KMC + 1)),
                  size = 2, shape = 21, fill = "#febf38", stroke = NA) +
  scale_x_continuous(limits = c(0, log10(max(df_TCGA$sums_TCGA_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  scale_y_continuous(limits = c(0, log10(max(df_TCGA$sums_TCGA_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    legend.position = "none"
  ) +
  scale_fill_gradientn(colours = rev(colors))

# Save TCGA scatterplot
ggsave(plot = p_all_tcga_filt, path = plots.dir,
       filename = paste0("scatterplot_", project, "_TCGA_FILT_all.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 6, bg = 'transparent')


# ---- CPTAC ----
# Mark k-mers that pass CPTAC filtering
df_CPTAC$cptac_filter <- ifelse(df_CPTAC$kmers %in% cptac_filter$kmers, "passed", "failed")

# Calculate density for scatterplot visualization
df_CPTAC$density_CPTAC <- get_density(
  log10(df_CPTAC$sums_CPTAC_Cancer_KMC + 1),
  log10(df_CPTAC$sums_CPTAC_Healthy_KMC + 1)
)

# Scatterplot: all CPTAC vs filtered highlighted
p_all_CPTAC_filt <- ggplot() +
  geom_point_rast(data = df_CPTAC,
                  mapping = aes(x = log10(sums_CPTAC_Healthy_KMC + 1),
                                y = log10(sums_CPTAC_Cancer_KMC + 1),
                                fill = log10(density_CPTAC)),
                  size = 2, shape = 21, stroke = NA) +
  geom_point_rast(data = cptac_filter,
                  mapping = aes(x = log10(sums_CPTAC_Healthy_KMC + 1),
                                y = log10(sums_CPTAC_Cancer_KMC + 1)),
                  size = 2, shape = 21, fill = "#febf38", stroke = NA) +
  scale_x_continuous(limits = c(0, log10(max(df_CPTAC$sums_CPTAC_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  scale_y_continuous(limits = c(0, log10(max(df_CPTAC$sums_CPTAC_Cancer_KMC) * 1.1)),
                     breaks = 0:9) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    legend.position = "none"
  ) +
  scale_fill_gradientn(colours = rev(colors))

# Save CPTAC scatterplot
ggsave(plot = p_all_CPTAC_filt, path = plots.dir,
       filename = paste0("scatterplot_", project, "_CPTAC_FILT_all.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 6, bg = 'transparent')

# ----------------------------
# INTERSECT FILTERED K-MERS
# ----------------------------
intersected <- intersect(cptac_filter$kmers, tcga_filter$kmers)
intersected_df <- df |> dplyr::filter(kmers %in% intersected)
# Create a list for Venn diagram
x <- list(
  CPTAC = cptac_filter$kmers,
  TCGA  = tcga_filter$kmers
)

# ----------------------------
# VENN DIAGRAM
# ----------------------------
ven <- ggvenn(data = x, fill_color = c("#40af75", "#1e6bb0"))
ven

# Save Venn diagram
ggsave(plot = ven, path = plots.dir, 
       filename = paste0("venn_", project, "_CPTAC_TCGA.png"),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot = ven, path = plots.dir, 
       filename = paste0("venn_", project, "_CPTAC_TCGA.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# SAVE INTERSECTED K-MERS
# ----------------------------
write(
  x = intersected,  
  file = paste0("~/Desktop/Cancer_paper/data/intersected_kmers_",project,".txt")
)

# ----------------------------
# ADDITIONAL FILTERING USING SRA DATA
# ----------------------------

# Set a cutoff for SRA counts to filter out common sequences
sra_cutoff <- 200

# Keep only intersected k-mers that are not abundant in SRA
cancer_specific_df <- intersected_df |>
  filter(sums_SRA < sra_cutoff)

# ----------------------------
# DENSITY PLOT (SRA COUNTS)
# ----------------------------
dens_plot_SRA <- ggplot() +
  geom_density(data = intersected_df,
               mapping = aes(x = log10(sums_SRA+1)),
               fill = "#dacdd9") +
  geom_vline(xintercept = log10(sra_cutoff),
             linetype = "dashed", linewidth = 1) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(r = 20)),
    axis.title.y = element_text(size = 20, margin = margin(r = 20)),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm")
  ) +
  xlab(label = "log10 SRA counts") +
  xlim(c(0, 1.5 * max(log10(intersected_df$sums_SRA)))) +
  ggtitle(label = paste0(length(cancer_specific_df$kmers),
                         " kmers have < ", sra_cutoff, " SRA counts"))

# Save density plots
ggsave(plot = dens_plot_SRA, path = plots.dir,
       filename = paste0("density_", project, "_SRA_counts.png"),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot = dens_plot_SRA, path = plots.dir,
       filename = paste0("density_", project, "_SRA_counts.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# VIOLIN PLOT (DISTRIBUTION OF COUNTS)
# ----------------------------
# Reshape into long format for plotting across datasets
tmp <- cancer_specific_df |>
  pivot_longer(cols = c(sums_CPTAC_Cancer,
                        sums_CPTAC_Healthy,
                        sums_TCGA_Cancer,
                        sums_TCGA_Healthy,
                        sums_SRA))

# Clean up dataset names for plotting
a <- tmp$name |>
  str_split_i(pattern = "sums_", i = 2) |>
  str_replace(pattern = "_", replacement = " ") |>
  str_replace(pattern = "Healthy", replacement = "Adjacent Normal")

tmp$name <- ifelse(a == "SRA", "Non-Cancer SRA", a)
tmp$name <- factor(tmp$name,
                   levels = c("CPTAC Cancer","TCGA Cancer",
                              "CPTAC Adjacent Normal","TCGA Adjacent Normal",
                              "Non-Cancer SRA"))

# Assign cancer vs non-cancer grouping
tmp$condition <- ifelse(tmp$name %in% c("CPTAC Cancer","TCGA Cancer"),
                        "cancer", "non-cancer")

# Plot distributions with violin + jitter
violin_counts <- ggplot() +
  geom_violin(data = tmp,
              mapping = aes(x = name, y = log10(value+1), fill = condition)) +
  geom_jitter(data = tmp,
              mapping = aes(x = name, y = log10(value+1)),
              size = 0.2, position = position_jitter(0.2)) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, margin = margin(r = 20)),
    axis.text = element_text(size = 18, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_line(linewidth = 1),
    plot.title = element_text(size = 18),
    axis.ticks.length = unit(.2, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size =18),
    axis.text.x = element_text(angle = 45, vjust = 0.7)
  ) +
  scale_fill_manual(values = c("#ce3701","#72cacf")) +
  ylab(label = "log10 counts") +
  scale_x_discrete(labels = function(x){sub("\\s", "\n", x)})

# Save violin plots
ggsave(plot = violin_counts, path = plots.dir,
       filename = paste0("violin_", project, "_candidate_counts.png"),
       dpi = 600, device = "png", width = 8, height = 6)
ggsave(plot = violin_counts, path = plots.dir,
       filename = paste0("violin_", project, "_candidate_counts.pdf"),
       dpi = 600, device = "pdf", width = 8, height = 6)

# ----------------------------
# SAVE FILTERED RESULTS
# ----------------------------
write_feather(cancer_specific_df,
              sink = paste0(analysis.dir, "/cancer_enriched_filtered_",
                            kmer_lenght,"mers_", project, ".arrow"))

write_delim(x = cancer_specific_df,
            file = paste0(analysis.dir, "/cancer_enriched_filtered_",
                          kmer_lenght,"mers_",project,".tsv"),
            delim = "\t")

# ----------------------------
# MERGE K-MERS INTO LONGER SEQUENCES
# ----------------------------
system(paste0("/Users/panagiotiskalogeropoulos/Desktop/dekupl-mergeTags/mergeTags -k 17 -m 8 -n ",
              analysis.dir, "/cancer_enriched_filtered_",kmer_lenght,"mers_",project,".tsv > ",
              analysis.dir,"/cancer_enriched_filtered_sequences_",project,".tsv"))

# Load merged candidate sequences
cancer_specific_sequences <- read.delim(
  paste0(analysis.dir,"/cancer_enriched_filtered_sequences_",project,".tsv")
)


