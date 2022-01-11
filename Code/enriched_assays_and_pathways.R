########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 26th July 2021
#
#
# Description:
#
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

# ## ToxPrints for the CATMoS chemicals
# catmos_toxprints <- read_tsv(
#   here::here(Input, "CATMoS_toxprints.tsv"),
#   col_types = cols())

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

# ## CTD chem-pathway enriched associations
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


# Create dataframes ready for chemotype enrichment ---------------------------     ---------------------------     ---------------------------


## Toxic chemicals
acute_toxicants <- acute_proc %>%
  mutate(toxic = if_else(nontoxic, 0, 1)) %>%
  select(CASRN, toxic)

## Very toxic chemicals
very_toxic <- acute_proc %>%
  select(CASRN, very_toxic) %>%
  mutate(very_toxic = if_else(very_toxic, 1, 0))

## Very toxic vs nontoxic chemicals
v_tox_vs_nontox <- acute_proc %>%
  filter(very_toxic == TRUE | nontoxic == TRUE) %>%
  mutate(very_toxic = if_else(very_toxic, 1, 0)) %>%
  select(CASRN, very_toxic)


# Enriched ToxCast assays - toxic chemicals ---------------------------     ---------------------------     ---------------------------


if (file.exists(here::here(Output, "toxic_enriched_assays.tsv"))){
  toxic_enriched_assays <- readr::read_tsv(
    here::here(Output, "toxic_enriched_assays.tsv"),
    col_types = cols()
  )
} else {

  ## ID ToxCast assays enriched for activity in toxic acute chemicals
  toxic_enriched_assays <- acute_toxicants %>%
    inner_join(assay_hitcalls, by = c("CASRN" = "casrn")) %>%
    select(-CASRN) %>%
    chemotype_enrichment(., endpoint_of_interest = toxic) %>%
    filter(p.value <= 0.05) %>%
    arrange(desc(OR), p.value)

  readr::write_tsv(
    toxic_enriched_assays,
    here::here(Output, "toxic_enriched_assays.tsv")
  )

}


# Enriched ToxCast assays - very toxic chemicals ---------------------------     ---------------------------     ---------------------------


if (file.exists(here::here(Output, "very_toxic_enriched_assays.tsv"))){
  very_toxic_enriched_assays <- readr::read_tsv(
    here::here(Output, "very_toxic_enriched_assays.tsv"),
    col_types = cols()
  )
} else {

  ## ID ToxCast assays enriched for activity in very toxic acute chemicals
  very_toxic_enriched_assays <- very_toxic %>%
    inner_join(assay_hitcalls, by = c("CASRN" = "casrn")) %>%
    select(-CASRN) %>%
    chemotype_enrichment(., endpoint_of_interest = very_toxic) %>%
    filter(p.value <= 0.05) %>%
    arrange(desc(OR), p.value)

  readr::write_tsv(
    very_toxic_enriched_assays,
    here::here(Output, "very_toxic_enriched_assays.tsv")
  )

}


# Enriched ToxCast assays - v.toxic vs nontoxic chemicals ---------------------------     ---------------------------     ---------------------------


if (file.exists(here::here(Output, "very_toxic_vs_nontoxic_enriched_assays.tsv"))){
  v_toxic_vs_nontox_enriched_assays <- readr::read_tsv(
    here::here(Output, "very_toxic_vs_nontoxic_enriched_assays.tsv"),
    col_types = cols()
  )
} else {

  ## ID ToxCast assays enriched for activity in very toxic acute chemicals
  v_toxic_vs_nontox_enriched_assays <- v_tox_vs_nontox %>%
    inner_join(assay_hitcalls, by = c("CASRN" = "casrn")) %>%
    select(-CASRN) %>%
    chemotype_enrichment(., endpoint_of_interest = very_toxic) %>%
    filter(p.value <= 0.05) %>%
    arrange(desc(OR), p.value)

  readr::write_tsv(
    v_toxic_vs_nontox_enriched_assays,
    here::here(Output, "very_toxic_vs_nontoxic_enriched_assays.tsv")
  )

}


# # Enriched CTD Pathways - toxic chemicals ---------------------------     ---------------------------     ---------------------------
#
#
# if (file.exists(here::here(Output, "toxic_enriched_pathways.tsv"))){
#   toxic_enriched_pathways <- readr::read_tsv(
#     here::here(Output, "toxic_enriched_pathways.tsv"),
#     col_types = cols()
#   )
# } else {
#
#   ## ID CTD Pathways enriched for activity in toxic acute chemicals
#   toxic_enriched_pathways <- acute_toxicants %>%
#     inner_join(adjusted_ctd_pathway_mat, by = c("CASRN" = "casrn")) %>%
#     select(-CASRN) %>%
#     chemotype_enrichment(., endpoint_of_interest = toxic) %>%
#     filter(p.value <= 0.05) %>%
#     arrange(desc(OR), p.value)
#
#   readr::write_tsv(
#     toxic_enriched_pathways,
#     here::here(Output, "toxic_enriched_pathways.tsv")
#   )
#
# }
#
#
# # Enriched CTD Pathways - very toxic chemicals ---------------------------     ---------------------------     ---------------------------
#
#
# if (file.exists(here::here(Output, "very_toxic_enriched_pathways.tsv"))){
#   very_toxic_enriched_pathways <- readr::read_tsv(
#     here::here(Output, "very_toxic_enriched_pathways.tsv"),
#     col_types = cols()
#   )
# } else {
#
#   ## ID CTD Pathways enriched for activity in very toxic acute chemicals
#   very_toxic_enriched_pathways <- very_toxic %>%
#     inner_join(adjusted_ctd_pathway_mat, by = c("CASRN" = "casrn")) %>%
#     select(-CASRN) %>%
#     chemotype_enrichment(., endpoint_of_interest = very_toxic) %>%
#     filter(p.value <= 0.05) %>%
#     arrange(desc(OR), p.value)
#
#   readr::write_tsv(
#     very_toxic_enriched_pathways,
#     here::here(Output, "very_toxic_enriched_pathways.tsv")
#   )
#
# }
#
#
# # Enriched CTD Pathways - v.toxic vs nontoxic chemicals ---------------------------     ---------------------------     ---------------------------
#
#
# if (file.exists(here::here(Output, "very_toxic_vs_nontoxic_enriched_pathways.tsv"))){
#   v_toxic_vs_nontox_enriched_pathways <- readr::read_tsv(
#     here::here(Output, "very_toxic_vs_nontoxic_enriched_pathways.tsv"),
#     col_types = cols()
#   )
# } else {
#
#   ## ID CTD Pathways enriched for activity in very toxic acute chemicals
#   v_toxic_vs_nontox_enriched_pathways <- v_tox_vs_nontox %>%
#     inner_join(adjusted_ctd_pathway_mat, by = c("CASRN" = "casrn")) %>%
#     select(-CASRN) %>%
#     chemotype_enrichment(., endpoint_of_interest = very_toxic) %>%
#     filter(p.value <= 0.05) %>%
#     arrange(desc(OR), p.value)
#
#   readr::write_tsv(
#     v_toxic_vs_nontox_enriched_pathways,
#     here::here(Output, "very_toxic_vs_nontoxic_enriched_pathways.tsv")
#   )
#
# }
#
#
