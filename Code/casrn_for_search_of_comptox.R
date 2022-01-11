########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 13th April 2020
#
#
# Description: Set ids up for searching through the CompTox Chemicals Dashboard
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
  tidyverse
)


Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


# Load data files ---------------------------     ---------------------------     ---------------------------


## CATMoS Processed Acute tox values
load(here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- as_tibble(fdat)


# Write files to search CompTox Dashboard ---------------------------     ---------------------------     ---------------------------


## Extract CATMoS CASRNs
acute_proc %>%
  select(CASRN) %>%
  write_tsv(here(Inter, "CATMoS_CAS.tsv"), col_names = FALSE)
