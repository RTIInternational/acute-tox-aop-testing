########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 14th June 2021
#
#
# Description: Use ToxPrint fingerprints to:
#             1. calculate Tanimoto distance (i.e. 1-similarity),
#             2. Cluster substances using "best" method
#             3. Write out substances in each cluster to separate folders
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
  readxl,
  scales,
  ggplot2,
  ggthemes,
  class,
  corrr,
  philentropy,
  cluster,
  factoextra,
  kableExtra,
  fingerprint,
  ggdendro,
  clValid,
  clv,
  rcdk,
  tidymodels,
  tidyverse
)


Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"

## Set seed for reproducibility
set.seed(2020)


# Check if necessary files exist before continuing ---------------------------     ---------------------------     ---------------------------


## Has the list of CASRNs been created to search for structure info from CompTox Dashboard?
if (!file.exists(here(Inter, "CATMoS_CAS.tsv"))) {
  source(here("Code/casrn_for_search_of_comptox.R"))

  stop("Please use 'IntermediateData/CATMoS_CAS.tsv' to retrieve structure information
       from https://comptox.epa.gov/dashboard/ and save it in the InputData folder")
}

## Has structure info been retrieved from the CompTox Dashboard?
if (!file.exists(list.files(here::here(Input), pattern = "CompToxChemicalsDashboard-Batch-Search.*tsv", full.names = TRUE))) {
  stop("Please use 'IntermediateData/CATMoS_CAS.tsv' to retrieve structure information
       from https://comptox.epa.gov/dashboard/ and save it in the InputData folder")
}

## Has the structure info from CompTox Dashboard been converted to SDF and
## has a ToxPrint fingerprint file been generated using the ChemoTyper app?
if (!file.exists(here(Inter, "CATMoS_SMILES.sdf")) | !file.exists(here(Input, "CATMoS_toxprints.tsv"))) {
  source(here("Code/smiles_to_sdf_for_toxprints.R"))
}


# Load data files ---------------------------     ---------------------------     ---------------------------


## CATMoS Processed Acute tox values
load(here::here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- tibble::as_tibble(fdat)

## SMILES info for CATMoS chemicals
catmos_smiles <- list.files(
  path = here::here(Input),
  pattern = "CompToxChemicalsDashboard-Batch-Search.*tsv",
  full.names = TRUE
) %>%
  purrr::set_names(.) %>%
  # Read each file and append to bottom of each file
  purrr::map_dfr(
    read_tsv,
    col_types = cols(),
    na = c("", "NA", "-", "NO_MATCH")
  ) %>%
  dplyr::filter(!is.na(QSAR_READY_SMILES)) %>%
  rename(CASRN = INPUT) %>%
  select(-FOUND_BY)

## ToxPrints of CATMoS chemicals
catmos_toxprints <- read_tsv(
  here::here(Input, "CATMoS_toxprints.tsv"),
  col_types = cols()
)


# Set parameters ---------------------------     ---------------------------     ---------------------------


## Which fingerprint currently being used?
fprint <- "Toxprint"

## What clustering algorithm to use?
clus_algo <- "ward"

## What cut height to use?
c_height <- 0.7


# Set up fingerprint folders ---------------------------     ---------------------------     ---------------------------
###   Set up folders for different fingerprints in    ###
###   IntermediateData and OutputFiles folder    ###


## If subfolder for fingerprint don't exist in
## IntermediateData create it
if (!dir.exists(here::here(Inter, fprint))) {
  # If fingerprint folder doesn't exist - create it
  dir.create(here::here(Inter, paste0(fprint)))
}

## If subfolder for fingerprint don't exist in
## OutputFiles create it
if (!dir.exists(here::here(Output, fprint))) {
  # If fingerprint folder doesn't exist - create it
  dir.create(here::here(Output, paste0(fprint)))
}


# Calculate Tanimoto distance ---------------------------     ---------------------------     ---------------------------
###   Calculate Tanimoto distances between all CATMoS chemicals   ###


## Tanimoto distance of all CATMoS chemicals with ToxPrints
if (file.exists(here::here(Output, fprint, "catmos_SMILES_Tani_sim.rds"))) {
  catmos_dist <- readr::read_rds(
    here::here(Output, fprint, "catmos_SMILES_Tani_sim.rds")
  )
} else {
  catmos_dist <- catmos_toxprints %>%
    tibble::column_to_rownames("CASRN") %>% # Make DTXSID col into rownames
    # Calculate Tanimoto similarity
    philentropy::distance(.,
                          method = "tanimoto",
                          use.row.names = TRUE
    )
  # Save as Rdata file
  saveRDS(
    catmos_dist,
    here::here(Output, fprint, "catmos_SMILES_Tani_sim.rds")
  )
}


# Cluster substances ---------------------------     ---------------------------     ---------------------------
###   Use Ward's method to cluster substances   ###


### Generate hierarchical clusters

## Generation of hierarchical cluster using Ward's takes ~16hrs on laptop
if (file.exists(here::here(Output, fprint, "catmos_ward_clus.rds"))) {
  catmos_hc_ward <- readr::read_rds(
    here::here(Output, fprint, "catmos_ward_clus.rds")
  )
} else {
  # Generate hierarchical clusters using Ward's method
  catmos_hc_ward <- cluster::agnes(as.dist(catmos_dist), method = clus_algo)
  # Save Rdata file
  saveRDS(
    catmos_hc_ward,
    here::here(Output, fprint, "catmos_ward_clus.rds")
  )
}

# dend <- ggdendro::dendro_data(catmos_hc_ward, type = "rectangle")
#
# ## Plot dendrogram
# ggplot(dend$segments) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_hline(yintercept = c_height, colour = "#4B9CD3", linetype = "dashed") +
#   geom_text(
#     data = dend$labels,
#     aes(x, y, label = label),
#     hjust = 1,
#     angle = 90
#   ) +
#   labs(title = paste(stringr::str_to_title(clus_algo), "Clustering")) +
#   scale_x_discrete(name = "Chemical") +
#   scale_y_continuous(
#     name = "Height",
#     limits = c(0, NA),
#     breaks = seq(0, 20, by = 2),
#     labels = scales::label_number()
#   ) +
#   ggthemes::theme_fivethirtyeight() +
#   theme(
#     axis.title = element_text(size = 10),
#     axis.text.x = element_blank(),
#     panel.grid.major.x = element_blank(),
#     legend.position = "none",
#     text = element_text(family = "Inconsolata")
#   )
#
# ggsave(
#   filename = here(Output, "dendrogram_of_chemicals.png"),
#   height = 20,
#   width = 20,
#   units = "cm",
#   device = "png"
# )


# Find cut height to use ---------------------------     ---------------------------     ---------------------------


## Create variable with CASRN in
num_clusters <- catmos_toxprints %>%
  dplyr::select(CASRN)

## Cut heights between 0.4 & 2
cuts <- seq(0.4, 2, by = 0.1)

## Loop through sequence of heights,
## cut cluster tree, and add to num_clusters variable
for (i in cuts) {
  col_name <- paste0("Height_", i)

  cut_clus <- cutree(as.hclust(catmos_hc_ward),h = i) %>%
    dplyr::as_tibble(rownames = "CASRN") %>%
    dplyr::rename(!!col_name := value)

  num_clusters <- dplyr::inner_join(num_clusters, cut_clus, by = "CASRN")
}

## Calculate:
## 1) number of chemicals per cluster per height cut,
## 2) number of single chemical clusters per height,
## 3) average chemicals per cluster per height, and
## 4) max chemicals in a cluster per height
cluster_stats <- num_clusters %>%
  tidyr::pivot_longer(cols = c(contains("Height"))) %>%
  dplyr::group_by(name, value) %>%
  dplyr::mutate(
    num_of_chems = n()
  ) %>%
  dplyr::distinct(name, num_of_chems) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(
    total_clusters = n_distinct(value),
    num_single_chem_clus = sum(num_of_chems == 1),
    num_pairs_chem_clus = sum(num_of_chems == 2),
    average = ave(num_of_chems),
    median = median(num_of_chems),
    max_chems = max(num_of_chems)
  )

## Show only info for cut heights
## not individual clusters
cluster_stats %>%
  dplyr::distinct(
    name,
    total_clusters,
    num_single_chem_clus,
    num_pairs_chem_clus,
    average,
    median,
    max_chems
  )


# Write out chemicals to their clusters ---------------------------     ---------------------------     ---------------------------
###   For each cut height write out tsv for each cluster    ###


## Cut heights between 0.7 & 1.4
cuts <- seq(0.7, 1.4, by = 0.1)

## If directory for cluster height of interest doesn't exist
## Work through list of cut heights and create separate folder for each
## cut height and a separate tsv for each cluster
if (!dir.exists(here::here(Output, fprint, paste0("Clusters_at_Height_", c_height)))) {

  for (i in cuts) {
    # If folder for cut height doesn't exist create one
    if (!dir.exists(here::here(Output, fprint, paste0("Clusters_at_Height_", i)))) {
      clus_height <- dir.create(
        here::here(Output, fprint, paste0("Clusters_at_Height_", i))
      )
    }

    # Cut tree at given height and join with SMILES info
    clusters <- cutree(as.hclust(catmos_hc_ward), h = i) %>%
      dplyr::as_tibble(rownames = "CASRN") %>%
      dplyr::rename(cluster_id = value) %>%
      dplyr::inner_join(catmos_smiles, by = "CASRN") %>%
      dplyr::inner_join(
        acute_proc %>%
          select(CASRN, very_toxic, nontoxic),
        by = "CASRN"
      )

    # Get IDs for each cluster
    group_ids <- dplyr::distinct(clusters, cluster_id) %>%
      purrr::as_vector()

    # For each cluster at given height create tibble of
    # SMILES, DTXSID, chem name, and CASRN & write out tsv
    clusters %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::arrange(cluster_id, DTXSID) %>%
      dplyr::select(SMILES, cluster_id, CASRN, DTXSID, PREFERRED_NAME, very_toxic, nontoxic) %>%
      dplyr::group_walk(~ readr::write_tsv(
        .x,
        here::here(Output, fprint, paste0("Clusters_at_Height_", i), paste0("Cluster_", .y$cluster_id, ".tsv")),
        col_names = FALSE
      ))
  }
}


# Write out cluster assignment w/CATMoS data ---------------------------     ---------------------------     ---------------------------


## Cut heights between 0.7 & 1
cuts <- seq(0.7, 1, by = 0.1)
cuts <- sort(cuts, decreasing = TRUE)

## Function to get cluster ID for different cut heights
clustering_assignments <- function(df, cut_height) {
  column_name <- paste0("cluster_id_", cut_height)

  clus_tmp <- cutree(as.hclust(df), h = cut_height) %>%
    dplyr::as_tibble(rownames = "CASRN") %>%
    dplyr::rename({{ column_name }} := value)
}

## Make copy of acute_proc
clus_assignments <- acute_proc

## For each cut height get cluster id for each chemicals
for (i in cuts) {

  clus_assignments <- clus_assignments %>%
    inner_join(clustering_assignments(catmos_hc_ward, i), by = "CASRN")

}

## Tibble containing CATMoS chemicals & the cluster each was assigned to
clus_assignments <- catmos_smiles %>%
  dplyr::inner_join(clus_assignments, by = "CASRN")


# Calculate Intra- & Inter-cluster statistics - all cut heights ---------------------------     ---------------------------     ---------------------------


## For each cut height calculate intra-cluster distances
for (i in cuts) {

  cluster_name <- paste0("cluster_id_", i)
  avg_col_name <- paste0("intracls_avg_", i)

  clusters <- clus_assignments %>% pull({{ cluster_name }})

  ## Calculate a number of intra- and inter-cluster stats
  intra_inter_clus_stats <- clv::cls.scatt.diss.mx(catmos_dist, clusters)

  ## Write out all stats as RData file
  write_rds(
    intra_inter_clus_stats,
    here(Output, fprint, paste0("intra_inter_clus_stats_cut_height_", i, ".rds"))
  )

  ## Get average intra-cluster distances
  avg_intraclus <- intra_inter_clus_stats[["intracls.average"]] %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "cluster_id", values_to = "intracls_avg") %>%
    mutate(
      cluster_id = str_remove_all(cluster_id, "c"),
      cluster_id = as.integer(cluster_id)
    ) %>%
    rename(
      {{ cluster_name }} := cluster_id,
      {{ avg_col_name }} := intracls_avg
    )


  clus_assignments <- clus_assignments %>%
    left_join(avg_intraclus, by = {{ cluster_name }}) %>%
    relocate({{ avg_col_name }}, .after = {{ cluster_name }})
}


## Write out cluster assignments with intra-cluster distances
if (!file.exists(here::here(Output, fprint, "Cluster_assignment.tsv"))) {
  write_tsv(
    clus_assignments,
    here::here(Output, fprint, "Cluster_assignment.tsv")
  )
}


# Calculate Inter-cluster statistics - cut height of interest ---------------------------     ---------------------------     ---------------------------


## Calculate intra- and inter-cluster stats for
## clusters at a specific cut-height
if (file.exists(here::here(Output, fprint, paste0("average_inter_clus_dist_cut_height_", c_height, ".rds")))){
  intra_inter_clus_stats <- read_rds(
    here::here(Output, fprint, paste0("average_inter_clus_dist_cut_height_", c_height, ".rds"))
  )
} else {

  cluster_name <- paste0("cluster_id_", c_height)
  clusters <- clus_assignments %>% pull({{ cluster_name }})

  ## Calculate intra- and inter-cluster stats
  intra_inter_clus_stats <- clv::cls.scatt.diss.mx(catmos_dist, clusters)

  ## Get just avg inter-cluster distances
  avg_interclus <- intra_inter_clus_stats[["intercls.average"]] %>%
    as_tibble(rownames = "cluster_id_1") %>%
    pivot_longer(cols = -cluster_id_1, names_to = "cluster_id_2", values_to = "intercls_avg") %>%
    mutate(
      cluster_id_1 = str_remove_all(cluster_id_1, "c"),
      cluster_id_1 = as.integer(cluster_id_1),
      cluster_id_2 = str_remove_all(cluster_id_2, "c"),
      cluster_id_2 = as.integer(cluster_id_2),
      intercls_avg = round(intercls_avg, 3)
    )

  ## Write avg inter-cluster distance as separate file
  write_tsv(
    avg_interclus,
    here(Output, fprint, paste0("average_inter_clus_dist_cut_height_", c_height, ".tsv"))
  )

}
