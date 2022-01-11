########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 28th Sept 2021
#
#
# Description:
#
# Notes:
#
#
# Potential Issues: None known
#########################################################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  janitor,
  patchwork,
  tidyverse
)


Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"

figs <- "Figures"


# Load data ---------------------------     ---------------------------     ---------------------------


## Load results from run_cluster_toxcast.R or run script to generate files
if(!file.exists(here(Output, "minimal_cluster_assays_summary_filtered_assays.tsv"))) {
  source("Code/run_cluster_toxcast.R")
} else {
  toxcast_chems <- read_tsv(here(Output, "toxcast_chems_filtered_assays.tsv"))

  chem_tox_activity_summary <- read_tsv(here(Output, "chem_tox_activity_summary_filtered_assays.tsv"))

  minimal_assays <- read_rds(here(Output, "minimal_assays_without_using_structural_clusters_filtered_assays.rds"))

  unique_mc_assays <- read_rds(here(Output, "minimal_assays_with_using_structural_clusters_filtered_assays.rds"))

  minimal_cluster_assays_summary <- read_tsv(here(Output, "minimal_cluster_assays_summary_filtered_assays.tsv"))
}

## CATMoS Processed Acute tox values
load(here::here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- tibble::as_tibble(fdat)


# Create dataframes identifying toxic and non-toxic chemicals ---------------------------     ---------------------------     ---------------------------


## Non-toxic chemicals
not_toxic <- acute_proc %>%
  filter(nontoxic == TRUE) %>%
  pull(CASRN)

## Toxic chemicals
acute_toxicants <- acute_proc %>%
  filter(nontoxic == FALSE) %>%
  pull(CASRN)

## Very toxic chemicals
very_toxic <- acute_proc %>%
  filter(very_toxic == TRUE) %>%
  pull(CASRN)


# Coverage of CATMoS chemicals - ToxCast ---------------------------     ---------------------------     ---------------------------
## Get numbers for pyramid plot (Figure 2A)

## Number of unique acutely toxic CATMoS chemicals tested in ToxCast
num_toxic_atwg_toxcast <- toxcast_chems %>%
  filter(nontoxic == FALSE) %>%
  pull(tested)
# 2076

## Number of unique acutely toxic CATMoS chemicals active in at least 1 ToxCast assay
num_toxic_atwg_act_1_assay <- toxcast_chems %>%
  filter(nontoxic == FALSE) %>%
  pull(active)
# 1886

## Number of unique acutely toxic CATMoS chemicals active below cytotoxicity in at least 1 ToxCast assay
num_toxic_atwg_act_bel_cyto <- toxcast_chems %>%
  filter(nontoxic == FALSE) %>%
  pull(cyto_hitcall)
# 1627


### NOT acutely toxic chemicals

## Number of unique not acutely toxic CATMoS chemicals tested in ToxCast
num_nontoxic_atwg_toxcast <- toxcast_chems %>%
  filter(nontoxic == TRUE) %>%
  pull(tested)
# 1917

## Number of unique not acutely toxic CATMoS chemicals active in at least 1 ToxCast assay
num_nontoxic_atwg_act_1_assay <- toxcast_chems %>%
  filter(nontoxic == TRUE) %>%
  pull(active)
# 1677

## Number of unique not acutely toxic CATMoS chemicals active below cytotoxicity in at least 1 ToxCast assay
num_nontoxic_atwg_act_bel_cyto <- toxcast_chems %>%
  filter(nontoxic == TRUE) %>%
  pull(cyto_hitcall)
# 1412


# Create pyramid (population) plots - within ToxCast ---------------------------     ---------------------------     ---------------------------
## Numbers are after removing cytoxicity and "Caution" assays

## Number of chemicals tested, active, & active below cytotox
chem_subset_toxcast <- tibble(
  chem_subset = factor(c(rep("atwg_total", 2),
                         rep("atwg_in_toxcast", 2),
                         rep("atwg_act_1_assay", 2),
                         rep("atwg_act_below_cyto", 2)),
                       levels = c("atwg_act_below_cyto", "atwg_act_1_assay", "atwg_in_toxcast", "atwg_total")),
  acute_toxic = factor(rep(c("nontoxic", "toxic"), 4), levels = c("nontoxic", "toxic")),
  value = c(length(not_toxic), length(acute_toxicants),
            num_nontoxic_atwg_toxcast, num_toxic_atwg_toxcast,
            num_nontoxic_atwg_act_1_assay, num_toxic_atwg_act_1_assay,
            num_nontoxic_atwg_act_bel_cyto, num_toxic_atwg_act_bel_cyto)
)

## Answer from "Eeeeed" here https://stackoverflow.com/questions/14680075/simpler-population-pyramid-in-ggplot2
## helped create this plot
toxcast_tested_pyramid_plot <- ggplot(chem_subset_toxcast, aes(
  x = chem_subset,
  y = ifelse(acute_toxic == "nontoxic", -value, value),
  fill = acute_toxic
)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = abs(value)),
    #size = 6,
    family = "Inconsolata",
    hjust = ifelse(chem_subset_toxcast$acute_toxic == "nontoxic", 1.1, -0.1)
  ) +
  scale_x_discrete(
    name = "Chemical Subset",
    labels = c(
      "atwg_act_below_cyto" = "Active below \ncytotoxicity",
      "atwg_act_1_assay" = "Active in at least \n 1 ToxCast assay",
      "atwg_in_toxcast" = "ATWG chemicals \nin ToxCast",
      "atwg_total" = "ATWG chemicals with \ntoxic or nontoxic designation"
    )
  ) +
  scale_y_continuous(
    limits = max(chem_subset_toxcast$value) * c(-1, 1) * 1.1,
    labels = abs
  ) +
  coord_flip() +
  scale_fill_manual(
    name = "ATWG Designation",
    values = c("#4B9CD3", "#bb0000"),
    breaks = c("nontoxic", "toxic")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = 14),
    legend.position = "bottom",
    aspect.ratio = 0.3
  )

ggsave(
  toxcast_tested_pyramid_plot,
  filename = here(Output, figs, "ATWG_tested_in_ToxCast_pyramid_plot.png"),
  width = 30,
  height = 15,
  units = "cm",
  device = "png"
)


# Number of chemicals identified by assay analysis workflow ---------------------------     ---------------------------     ---------------------------
## Get numbers for pyramid plot (Figure 2B)


## Number of unique acutely toxic CATMoS chemicals tested identified
## in analysis workflow without using structural clustering
num_toxic_all_assays <- chem_tox_activity_summary %>%
  filter(toxic == "Yes") %>%
  pull(all_assay_positives)

## Number of unique acutely toxic CATMoS chemicals tested identified
## in analysis workflow using structural clustering
num_toxic_with_struc_clus <- chem_tox_activity_summary %>%
  filter(toxic == "Yes") %>%
  pull(cluster_assay_positives)

## Number of unique not acutely toxic CATMoS chemicals tested identified
## in analysis workflow without using structural clustering
num_nontoxic_all_assays <- chem_tox_activity_summary %>%
  filter(toxic == "No") %>%
  pull(all_assay_positives)

## Number of unique not acutely toxic CATMoS chemicals tested identified
## in analysis workflow using structural clustering
num_nontoxic_with_struc_clus <- chem_tox_activity_summary %>%
  filter(toxic == "No") %>%
  pull(cluster_assay_positives)


# Create pyramid (population) plots - Minimal assays ---------------------------     ---------------------------     ---------------------------


## Only looking at minimal assays
chem_subset_min_assays <- tibble(
  chem_subset = factor(c(rep("atwg_all_assays", 2),
                         rep("atwg_w_struc_specific_assay", 2)),
                       levels = c("atwg_w_struc_specific_assay", "atwg_all_assays")),
  acute_toxic = factor(rep(c("nontoxic", "toxic"), 2), levels = c("nontoxic", "toxic")),
  value = c(num_nontoxic_all_assays, num_toxic_all_assays,
            num_nontoxic_with_struc_clus, num_toxic_with_struc_clus)
)

## Answer from "Eeeeed" here https://stackoverflow.com/questions/14680075/simpler-population-pyramid-in-ggplot2
## helped create this plot
min_assay_pyramid_plot <- ggplot(chem_subset_min_assays, aes(
  x = chem_subset,
  y = ifelse(acute_toxic == "nontoxic", -value, value),
  fill = acute_toxic
)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = abs(value)),
    #size = 6,
    family = "Inconsolata",
    hjust = ifelse(chem_subset_min_assays$acute_toxic == "nontoxic", 1.1, -0.1)
  ) +
  scale_x_discrete(
    name = "Analysis Workflow",
    labels = c(
      "atwg_w_struc_specific_assay" = "With Structural Clustering",
      "atwg_all_assays" = "Without Structural Clustering"
    )
  ) +
  scale_y_continuous(
    limits = max(chem_subset_min_assays$value) * c(-1, 1) * 1.1,
    labels = abs
  ) +
  coord_flip() +
  scale_fill_manual(
    name = "ATWG Designation",
    values = c("#4B9CD3", "#bb0000"),
    breaks = c("nontoxic", "toxic")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = 14),
    legend.position = "bottom",
    aspect.ratio = 0.2
  )

ggsave(
  min_assay_pyramid_plot,
  filename = here(Output, figs, "ATWG_ToxCast_minimal_assay_pyramid_plot.png"),
  width = 30,
  height = 15,
  units = "cm",
  device = "png"
)

## Combine pyramid plots & save figure
toxcast_tested_pyramid_plot / min_assay_pyramid_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2, tag_level = "keep", guides = "collect") &
  theme(legend.position='bottom')

ggsave(
  filename = here(Output, figs, "Combined_ATWG_pyramid_plots.png"),
  width = 30,
  height = 32,
  units = "cm",
  device = "png"
)


# Bar plots for minimal assay count & FDR ---------------------------     ---------------------------     ---------------------------


## Calculate False Discovery Rates
## Without structural clustering
all_assay_fdr <- round(num_nontoxic_all_assays/(num_nontoxic_all_assays + num_toxic_all_assays), digits = 2) * 100

## With structural clustering
struc_clus_fdr <- round(num_nontoxic_with_struc_clus/(num_nontoxic_with_struc_clus + num_toxic_with_struc_clus), digits = 2) * 100

## Create dataframe for use in generating
assay_count_fdr <- tibble(
  chem_subset = factor(c(rep("without_struc_clus", 3),
                         rep("with_struc_clus", 3)),
                       levels = c("without_struc_clus", "with_struc_clus")),
  measure = factor(rep(c("assay_count", "fdr", "assay_cover_single_chem"), 2),
                   levels = c("assay_count", "fdr", "assay_cover_single_chem")),
  # Assay Count, FDR, Num of assays required to cover a single chemical
  value = c(length(minimal_assays), all_assay_fdr, length(minimal_assays),
            length(unique_mc_assays), struc_clus_fdr, max(minimal_cluster_assays_summary$n))
)

## Bar plot of assay counts from minimal assay sets
assay_count_plot <- assay_count_fdr %>%
  filter(measure == "assay_count") %>%
  ggplot(., aes(x = chem_subset, y = value)) +
  geom_col(fill = "#002147") +
  geom_text(aes(label = abs(value)),
            #size = 6,
            family = "Inconsolata",
            vjust = -1) +
  ggtitle("Minimum Number of Assays Required \nto Capture All Toxic Chemicals") +
  scale_x_discrete(
    labels = c(
      "without_struc_clus" = "Without Structural \nClustering",
      "with_struc_clus" = "With Structural \nClustering"
    )
  ) +
  scale_y_continuous(
    name = "Assay Count: All Chemicals",
    breaks = scales::pretty_breaks(n = 7),
    limits = c(0, 350)
  ) +
  ggthemes::theme_fivethirtyeight() +
  theme(
    axis.text.x = element_text(),
    axis.title.y = element_text(),
    text = element_text(size = 14),
    legend.position = "bottom"
  )

## Bar plot of FDR from minimal assay sets
fdr_plot <- assay_count_fdr %>%
  filter(measure == "fdr") %>%
  ggplot(., aes(x = chem_subset, y = value)) +
  geom_col(fill = "#002147") +
  geom_text(aes(label = abs(value)),
            #size = 6,
            family = "Inconsolata",
            vjust = -1) +
  ggtitle("False Discovery Rate When \nUsing Minimal Assay Set") +
  scale_x_discrete(
    labels = c(
      "without_struc_clus" = "Without Structural \nClustering",
      "with_struc_clus" = "With Structural \nClustering"
    )
  ) +
  scale_y_continuous(
    name = "False Discovery Rate (%)",
    breaks = scales::pretty_breaks(n = 7),
    limits = c(0, 50)
  ) +
  ggthemes::theme_fivethirtyeight() +
  theme(
    axis.text.x = element_text(),
    axis.title.y = element_text(),
    text = element_text(size = 14),
    legend.position = "bottom"
  )

## Combine Assay Count and FDR into plot with A and B tags
# assay_count_plot + fdr_plot +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(ncol = 2, tag_level = "keep", guides = "collect") &
#   theme(legend.position='bottom')
#
# ggsave(
#   filename = here(Output, figs, "Combined_assay_count_fdr_bar_plots.png"),
#   width = 47,
#   height = 30,
#   units = "cm",
#   device = "png"
# )


# Assays needed to cover single chemical ---------------------------     ---------------------------     ---------------------------


## Limit dataframe to assay_cover_single_chem measure
cover_single_chem <- assay_count_fdr %>%
  filter(measure == "assay_cover_single_chem") %>%
  add_column(graph_label = c("177", "< 5"))

## Create bar chart
cover_single_chem_plot <- ggplot(cover_single_chem, aes(x = chem_subset, y = value)) +
  geom_col(fill = "#002147") +
  geom_text(aes(label = graph_label),
            family = "Inconsolata",
            vjust = -1) +
  ggtitle("Number of Assays Needed to \nCover a Single Chemical") +
  scale_x_discrete(
    labels = c(
      "without_struc_clus" = "Without Structural \nClustering",
      "with_struc_clus" = "With Structural \nClustering"
    )
  ) +
  scale_y_continuous(
    name = "Assay Count",
    breaks = scales::pretty_breaks(n = 7),
    limits = c(0, max(cover_single_chem$value) + 15)
  ) +
  ggthemes::theme_fivethirtyeight() +
  theme(
    axis.text.x = element_text(),
    axis.title.y = element_text(),
    text = element_text(size = 14),
    legend.position = "bottom"
  )

ggsave(
  cover_single_chem_plot,
  filename = here(Output, figs, "assays_needed_to_cover_single_chemical.png"),
  width = 15,
  height = 15,
  units = "cm",
  device = "png"
)


# Cluster count vs assay count per cluster ---------------------------     ---------------------------     ---------------------------


## Use info from 'minimal_cluster_assays_summary' variable
assay_count_per_clus <- tibble(
  assay_count = factor(minimal_cluster_assays_summary$n, levels = c("1", "2", "3", "4")),
  cluster_count = minimal_cluster_assays_summary$nn
)


assay_count_per_clus_plot <- ggplot(assay_count_per_clus, aes(x = assay_count, y = cluster_count)) +
  geom_col(fill = "#002147") +
  geom_text(aes(label = cluster_count),
            family = "Inconsolata",
            vjust = -1) +
  ggtitle("Cluster Count vs Assay Count \nper Cluster") +
  scale_x_discrete(
    name = "Assay Count per Cluster"
  ) +
  scale_y_continuous(
    name = "Cluster Count",
    breaks = scales::pretty_breaks(n = 7),
    limits = c(0, max(assay_count_per_clus$cluster_count) + 15)
  ) +
  ggthemes::theme_fivethirtyeight() +
  theme(
    axis.text.x = element_text(),
    axis.title = element_text(),
    text = element_text(size = 14),
    legend.position = "bottom"
  )

ggsave(
  assay_count_per_clus_plot,
  filename = here(Output, figs, "cluster_count_vs_assay_count_per_cluster.png"),
  width = 25,
  height = 15,
  units = "cm",
  device = "png"
)

## Combine assays needed to assess single chemical and
## assay count vs cluster count plots
cover_single_chem_plot + assay_count_per_clus_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, widths = c(1, 1.75), tag_level = "keep", guides = "collect") &
  theme(legend.position='bottom')

ggsave(
  filename = here(Output, figs, "Combined_cover_single_chem_assay_count_per_clus_plots.png"),
  width = 40,
  height = 30,
  units = "cm",
  device = "png"
)


# Combine figures into 1 for Fig. 3 A-D ---------------------------     ---------------------------     ---------------------------


## Combine assays needed to assess single chemical and
## assay count vs cluster count plots
assay_count_plot + fdr_plot +
  cover_single_chem_plot + assay_count_per_clus_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, widths = c(1, 1), tag_level = "keep", guides = "collect") &
  theme(title = element_text(size = 18),
        text = element_text(size = 16))

ggsave(
  filename = here(Output, figs, "Combined_quad_plot_for_Figure_3.png"),
  width = 40,
  height = 30,
  units = "cm",
  device = "png"
)

