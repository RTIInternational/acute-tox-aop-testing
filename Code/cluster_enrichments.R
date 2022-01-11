########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 14th June 2021
#
#
# Description: Use ToxPrint fingerprints to calculate: :
#             1. Chemotypes enriched in Chemical clusters,
#             2. ToxCast assays enriched in chemical clusters, and
#             3. CTD Pathways enriched in chemical clusters
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
  svDialogs,
  tidymodels,
  tidyverse
)
source(here::here("Code/chemotype_enrich_func.R"))
#source(here::here("Code/ctd_pathways_phenos_and_diseases_to_remove.R"))


Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


# Load data files ---------------------------     ---------------------------     ---------------------------


## CATMoS Processed Acute tox values
load(here::here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- tibble::as_tibble(fdat)

## ToxPrints for the CATMoS chemicals
catmos_toxprints <- read_tsv(
  here::here(Input, "CATMoS_toxprints.tsv"),
  col_types = cols())

## ToxCast hitcalls
## Adjusted using the lower bound of the cytotoxicity point
assay_hitcalls <- readr::read_tsv(
  here::here(Inter, "lower_cyto_pt_adjusted_hitcalls_all_assays_invitrodb_v3-3.tsv"),
  col_types = cols()
) %>%
  tidyr::pivot_wider(
    id_cols = "casrn",
    names_from = "assay_name",
    values_from = "cyto_hitcall"
  )

## Cluster assignments from hierarchical clustering
clus_assignments <- readr::read_tsv(
  here::here(Output, "Toxprint/Cluster_assignment.tsv"),
  col_types = cols()
)


## CTD chem-pathway enriched associations
# adjusted_ctd_pathway_mat <- read_tsv(
#   here(Input, "CTD/CTD_chem_pathways_enriched.tsv"),
#   skip = 27
# ) %>%
#   rename(
#     ChemicalName = `# ChemicalName`,
#     CASRN = CasRN
#   ) %>%
#   drop_na(CASRN) %>%
#   inner_join(
#     acute_proc %>%
#       select(CASRN),
#     by = "CASRN"
#   ) %>%
#   distinct(CASRN, PathwayID, .keep_all = TRUE) %>%
#   filter(
#     str_detect(PathwayName, paste0(path_pheno_to_remove, collapse = "|"), negate = TRUE)
#   ) %>%
#   distinct(CASRN, PathwayName, .keep_all = TRUE) %>%
#   mutate(effect = str_c("Pathway_", PathwayName),
#          value = if_else(CorrectedPValue <= mean(.$CorrectedPValue) + (1 * sd(.$CorrectedPValue)), 1, 0)) %>%
#   complete(CASRN, effect, fill = list(value = 0)) %>%
#   select(CASRN, effect, value) %>%
#   pivot_wider(
#     id_cols = CASRN,
#     names_from = effect,
#     values_from = value
#   ) %>%
#   rename(casrn = CASRN)

## CASRN for acute toxicants
toxic_casrn <- acute_proc %>%
  filter(nontoxic == FALSE) %>%
  pull(CASRN)


# Amend clus_assignments based on cut height of interest ---------------------------     ---------------------------     ---------------------------


#cut_height <- dlgInput("Enter a cut height between 0.7-1", Sys.info()[""])$res
cut_height <- 0.7
col_name_of_interest <- paste0("cluster_id_", cut_height)


# Chemotypes enriched within clusters ---------------------------     ---------------------------     ---------------------------
###   Using full ToxPrint name    ###


## If enrichments dataframe already exists load it
## else run code to generate enrichments and save it
# if (file.exists(here::here(Output, paste0("toxprints_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")))) {
#   cluster_toxprint_enrichments <- readr::read_tsv(
#     here::here(Output, paste0("toxprints_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")),
#     col_types = cols()
#   )
# } else {
# 
#   ## Link cluster assignments to ToxPrints
#   clusters <- clus_assignments %>%
#     dplyr::select(CASRN, clus_num = {{ col_name_of_interest }}) %>%
#     dplyr::inner_join(catmos_toxprints, by = "CASRN")
# 
#   ## Create vector of cluster numbers
#   ## for chemotype enrichment
#   cluster_ids <- clusters %>%
#     dplyr::distinct(clus_num) %>%
#     dplyr::pull()
# 
#   ## Set the max length of progress bar
#   ## to be the length of the cluster_ids vector
#   pb <- txtProgressBar(
#     min = 0,
#     max = length(cluster_ids),
#     initial = 0,
#     style = 3
#   )
# 
#   ## Create blank tibbles
#   cluster_toxprint_enrichments <- tibble::tibble()
#   tmp_enrich <- tibble::tibble()
# 
#   ## ID ToxPrints enriched within each cluster
#   for (i in seq_along(cluster_ids)) {
#     tmp_enrich <- clusters %>%
#       # If clus_num == current cluster ID set value
#       # in endpoint col as 1, otherwise set as 0
#       dplyr::mutate(
#         endpoint = if_else(clus_num == cluster_ids[[i]], 1, 0),
#         .after = clus_num
#       ) %>%
#       # Remove CASRN and clus_num cols
#       dplyr::select(-CASRN, -clus_num) %>%
#       # Perform chemotype enrichment
#       chemotype_enrichment(., endpoint_of_interest = endpoint) %>%
#       # Add current cluster number to tibble
#       dplyr::mutate(
#         cluster_id = cluster_ids[[i]],
#         .before = chemotype
#       ) %>%
#       # Keep only chemotype enrichments w/p-value <=0.05
#       dplyr::filter(p.value <= 0.05)
# 
#     # Add enrichments for current cluster to previous enrichments
#     cluster_toxprint_enrichments <- dplyr::bind_rows(cluster_toxprint_enrichments,
#                                                      tmp_enrich)
# 
#     # Make progress bar appear in Console
#     setTxtProgressBar(pb, i)
#   }
# 
#   write_tsv(
#     cluster_toxprint_enrichments,
#     here::here(Output, paste0("toxprints_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv"))
#   )
# }


# ToxCast assays enriched within clusters ---------------------------     ---------------------------     ---------------------------
###   Using adjusted hit calls based upon cytotoxicity   ###


## Extract hitcall data for CATMoS toxic chems in ToxCast assays
toxic_catmos_hitcalls <- assay_hitcalls %>%
  filter(casrn %in% toxic_casrn)

## If enrichments dataframe already exists load it
## else run code to generate enrichments and save it
if (file.exists(here::here(Output, paste0("toxcast_assays_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")))) {
  cluster_toxcast_enrichments <- readr::read_tsv(
    here::here(Output, paste0("toxcast_assays_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")),
    col_types = cols()
  )
} else {

  ## Link cluster info to acute tox-related ToxCast assays
  clusters <- clus_assignments %>%
    dplyr::select(CASRN, clus_num = {{ col_name_of_interest }}) %>%
    dplyr::inner_join(toxic_catmos_hitcalls, by = c("CASRN" = "casrn"))

  ## Create vector of cluster numbers
  ## for chemotype enrichment
  cluster_ids <- clusters %>%
    dplyr::distinct(clus_num) %>%
    dplyr::pull()

  ## Set the max length of progress bar
  ## to be the length of the cluster_ids vector
  pb <- txtProgressBar(
    min = 0,
    max = length(cluster_ids),
    initial = 0,
    style = 3
  )

  ## Create blank tibbles
  cluster_toxcast_enrichments <- tibble::tibble()
  tmp_enrich <- tibble::tibble()

  ## ID ToxCast assays enriched within each cluster
  for (i in seq_along(cluster_ids)) {
    tmp_enrich <- clusters %>%
      # If clus_num == current cluster ID set value
      # in endpoint col as 1, otherwise set as 0
      dplyr::mutate(
        endpoint = if_else(clus_num == cluster_ids[[i]], 1, 0),
        .after = clus_num
      ) %>%
      # Remove CASRN and clus_num cols
      dplyr::select(-CASRN, -clus_num) %>%
      # Perform chemotype enrichment
      chemotype_enrichment(., endpoint_of_interest = endpoint) %>%
      # Add current cluster number to tibble
      dplyr::mutate(
        cluster_id = cluster_ids[[i]],
        .before = chemotype
      ) %>%
      # Keep only toxcast assays w/p-value <=0.05
      dplyr::filter(p.value <= 0.05)

    # Add enrichments for current cluster to previous enrichments
    cluster_toxcast_enrichments <- dplyr::bind_rows(cluster_toxcast_enrichments,
                                                    tmp_enrich)

    # Make progress bar appear in Console
    setTxtProgressBar(pb, i)
  }

  write_tsv(
    cluster_toxcast_enrichments,
    here::here(Output, paste0("toxcast_assays_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv"))
  )
}


# CTD Pathways enriched within clusters ---------------------------     ---------------------------     ---------------------------
###   Using adjusted chem-pathway designations  ###


# ## Extract active Pathway data for CATMoS toxic chems
# toxic_catmos_pathways <- adjusted_ctd_pathway_mat %>%
#   filter(casrn %in% toxic_casrn)

# ## If enrichments dataframe already exists load it
# ## else run code to generate enrichments and save it
# if (file.exists(here::here(Output, paste0("CTD_pathways_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")))) {
#   cluster_pathway_enrichments <- readr::read_tsv(
#     here::here(Output, paste0("CTD_pathways_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv")),
#     col_types = cols()
#   )
# } else {
#
#   ## Link cluster info to CTD pathway
#   clusters <- clus_assignments %>%
#     dplyr::select(CASRN, clus_num = {{ col_name_of_interest }}) %>%
#     dplyr::inner_join(adjusted_ctd_pathway_mat, by = c("CASRN" = "casrn"))
#
#   ## Create vector of cluster numbers
#   ## for chemotype enrichment
#   cluster_ids <- clusters %>%
#     dplyr::distinct(clus_num) %>%
#     dplyr::pull()
#
#   ## Set the max length of progress bar
#   ## to be the length of the cluster_ids vector
#   pb <- txtProgressBar(
#     min = 0,
#     max = length(cluster_ids),
#     initial = 0,
#     style = 3
#   )
#
#   ## Create blank tibbles
#   cluster_pathway_enrichments <- tibble::tibble()
#   tmp_enrich <- tibble::tibble()
#
#   ## ID CTD Pathways enriched within each cluster
#   for (i in seq_along(cluster_ids)) {
#     tmp_enrich <- clusters %>%
#       # If clus_num == current cluster ID set value
#       # in endpoint col as 1, otherwise set as 0
#       dplyr::mutate(
#         endpoint = if_else(clus_num == cluster_ids[[i]], 1, 0),
#         .after = clus_num
#       ) %>%
#       # Remove CASRN and clus_num cols
#       dplyr::select(-CASRN, -clus_num) %>%
#       # Perform chemotype enrichment
#       chemotype_enrichment(., endpoint_of_interest = endpoint) %>%
#       # Add current cluster number to tibble
#       dplyr::mutate(
#         cluster_id = cluster_ids[[i]],
#         .before = chemotype
#       ) %>%
#       # Keep only pathway enrichments w/p-value <=0.05
#       dplyr::filter(p.value <= 0.05)
#
#     # Add enrichments for current cluster to previous enrichments
#     cluster_pathway_enrichments <- dplyr::bind_rows(cluster_pathway_enrichments,
#                                                     tmp_enrich)
#
#     # Make progress bar appear in Console
#     setTxtProgressBar(pb, i)
#   }
#
#   write_tsv(
#     cluster_pathway_enrichments,
#     here::here(Output, paste0("CTD_pathways_enriched_in_chemical_clusters_cut_height_", cut_height, ".tsv"))
#   )
# }


