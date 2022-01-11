########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0
#
#
# Description: Creates SDF using the SMILES from search of CompTox Chemicals Dashboard.
#             The SDF can be run in the ChemoTyper to generate ToxPrints fingerprints
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
  rcdk,
  tidyverse
)

Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


# Load data files ---------------------------     ---------------------------     ---------------------------


## QSAR-ready SMILES of CATMoS chemicals
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

## CATMoS Processed Acute tox values
load(here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- as_tibble(fdat)


# Write SDF to get Toxprints ---------------------------     ---------------------------     ---------------------------


## Parse SMILES so can create SDF for CATMoS chemicals
catmos_for_sdf <- catmos_smiles %>%
  dplyr::select(CASRN, QSAR_READY_SMILES) %>%
  tibble::deframe() %>%
  rcdk::parse.smiles()

## Create vector of CASRNs to be used as labels
catmos_casrn <- catmos_smiles %>%
  dplyr::distinct(CASRN) %>%
  purrr::as_vector()

## Assign CASRN as title of parsed SMILES
for (i in catmos_casrn) {
  rcdk::set.title(catmos_for_sdf[[i]], title = paste(i))
}

## Write structures out as SDF
rcdk::write.molecules(catmos_for_sdf, here::here(Inter, "CATMoS_SMILES.sdf"))


# Add chem name to CATMoS info ---------------------------     ---------------------------     ---------------------------


catmos_chem_info <- acute_proc %>%
  left_join(catmos_smiles, by = "CASRN") %>%
  select(DTXSID, PREFERRED_NAME, CASRN, very_toxic,
         nontoxic, LD50_mgkg, EPA_category, GHS_category)

write_tsv(catmos_chem_info, here(Inter, "full_catmos_with_chem_names.tsv"))


# Rename M_NAME col from ChemoTyper ---------------------------     ---------------------------     ---------------------------


#####                        #####                          #####
### CAN ONLY BE RUN AFTER TOXPRINTS CREATED USING CHEMOTYPER  ###
#####                        #####                          #####

# If it runs on a file that has been previously processed, you'll get the following error:
#  Error: Can't rename columns that don't exist.
#  x Column `M_NAME` doesn't exist.
#  Run `rlang::last_error()` to see where the error occurred.
# This is OK because the column has already been renamed
if (!file.exists(here(Input, "CATMoS_toxprints.tsv"))) {
  stop("Use 'IntermediateData/CATMoS_SMILES.sdf' and ChemoTyper app with ToxPrints v2.0_r711 to create ToxPrint
       fingerprint for CATMOS chemicals and save file in InputData folder as 'CATMOS_toxprints.tsv'")
} else {
  read_tsv(here::here(Input, "CATMoS_toxprints.tsv"), col_types = cols()) %>%
    rename(CASRN = M_NAME) %>%
    write_tsv(here::here(Input, "CATMoS_toxprints.tsv"))
}

