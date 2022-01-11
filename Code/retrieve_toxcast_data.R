########################################################################################
#
# Author: Mark D. Nelms, Ph.D., mnelms@rti.org
#
# Version: 1.0 16th June 2021
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
  RMariaDB, #Needed for queries not supported by tcpl
  data.table,
  tidyverse,
  configr
)

if (!require("tcpl")) {
  devtools::install_github("USEPA/CompTox-ToxCast-tcpl", ref = "dots_for_tcplConf")
  library(tcpl)
} else {
  library(tcpl)
}


Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


# Load data ---------------------------     ---------------------------     ---------------------------


## CATMoS Processed Acute tox values
load(here(Input, "AcuteTox FullDataset Processed 171130.RData"))

## CATMoS data as a tibble
acute_proc <- as_tibble(fdat)


# Connect to invitroDB v3.3 ---------------------------     ---------------------------     ---------------------------

##Read database configuration from yaml file
if (is.yaml.file(here("dbconfig.yml"))){
  dbconfig <- read.config(file = here("dbconfig.yml"))
}else{
  stop("Unable to read database configuration.")
}

## Connect to DB
tcplConf(drvr = "MySQL",
         user = dbconfig[['user']],
         pass = dbconfig[['password']],
         host = dbconfig[['host']],
         db = dbconfig[['database']],
         port = dbconfig[['dbport']])


# Retrieve chemical information ---------------------------     ---------------------------     ---------------------------


if (!file.exists(here::here(Inter, "chemical_information_invitrodb_v3-3.tsv"))) {
  chem_info <- tcplLoadChem(include.spid = FALSE)

  write_tsv(chem_info, here(Inter, "chemical_information_invitrodb_v3-3.tsv"))
} else {
  chem_info <- read_tsv(
    here(Inter, "chemical_information_invitrodb_v3-3.tsv"),
    na = c("", "NA", "-"),
    col_types = cols())
}


# Retrieve cytotoxicity information ---------------------------     ---------------------------     ---------------------------


if (file.exists(here::here(Inter, "cyto_point_invitrodb_v3-3.tsv"))){
  cyto_pt <- read_tsv(
    here(Inter, "cyto_point_invitrodb_v3-3.tsv"),
    col_types = cols()
  ) %>%
    select(casn, cyto_pt_um, lower_bnd_um)
} else {
  ## invitroDB v3.3 Cytotoxicity Info
  cyto_pt <- tcplCytoPt()

  write_tsv(cyto_pt, here(Inter, "cyto_point_invitrodb_v3-3.tsv"))
}


# Create chemical by assay matrices ---------------------------     ---------------------------     ---------------------------


if (file.exists(here(Inter, "ac50_matrix_all_assays_invitrodb_v3-3.tsv"))) {
  ac50_mat <- read_tsv(
    here(Inter, "ac50_matrix_all_assays_invitrodb_v3-3.tsv"),
    col_types = cols(casrn = col_character(), .default = col_double())
  )

  tested_mat <- read_tsv(
    here(Inter, "tested_matrix_all_assays_invitrodb_v3-3.tsv"),
    col_types = cols()
  )
} else {
  ## When more than one sample is included for a chemical/assay pair, tcplVarMat
  ## aggregates multiple samples to a chemical level call utilizing tcplSubsetChid.
  mats <- tcplVarMat(
    chid = NULL,
    aeid = NULL,
    add.vars = NULL,
    row.id = "casn", # options: casn, code, chnm, dsstox_substance_id
    flag = TRUE,
    cyto.pars = list(),
    include.na.chid = FALSE,
    odir = NULL,
    file.prefix = NULL
  )

  ## Create ac50 matrix
  ac50_mat <- as_tibble(mats$ac50, rownames = "casrn")

  write_tsv(ac50_mat, here(Inter, "ac50_matrix_all_assays_invitrodb_v3-3.tsv"))

  ## Create tested matrix
  tested_mat <- as_tibble(mats$tested, rownames = "casrn")

  write_tsv(tested_mat, here(Inter, "tested_matrix_all_assays_invitrodb_v3-3.tsv"))

  ## Create hitcall matrix
  hitc <- as_tibble(mats$hitc, rownames = "casrn")

  write_tsv(hitc, here(Inter, "hitcall_matrix_all_assays_invitrodb_v3-3.tsv"))
}


# Adjust hitcalls using cytotox data ---------------------------     ---------------------------     ---------------------------

if (file.exists(here(Inter, "lower_cyto_pt_adjusted_hitcalls_all_assays_invitrodb_v3-3.tsv"))) {
  ac50_mat_long <- read_tsv(
    here(Inter, "toxcast_activity_all_assays_invitrodb_v3-3.tsv"),
    col_types = cols(
      casrn = col_character(),
      assay_name = col_character(),
      .default = col_double()
      )
  )

  lwr_bnd_hitc <- read_tsv(
    here(Inter, "lower_cyto_pt_adjusted_hitcalls_all_assays_invitrodb_v3-3.tsv"),
    col_types = cols(
      casrn = col_character(),
      assay_name = col_character(),
      .default = col_double()
      )
  )
} else {
  ## Convert tested matrix to long form
  tested_mat_long <- tested_mat %>%
    pivot_longer(cols = -casrn, names_to = "assay_name", values_to = "tested")

  ac50_mat_long <- ac50_mat %>%
    #filter(casrn %in% acute_proc$CASRN) %>%
    pivot_longer(
      cols = -casrn,
      names_to = "assay_name",
      values_to = "assay_ac50"
    ) %>%
    mutate(assay_ac50 = replace_na(assay_ac50, replace = 1000000)) %>%
    inner_join(tested_mat_long, by = c("casrn", "assay_name"))

  ## Use lower cytotox point to create hitcall matrix adjusted
  ## for cytotoxicity
  lwr_bnd_hitc <- ac50_mat_long %>%
    inner_join(cyto_pt, by = c("casrn" = "casn")) %>%
    mutate(cyto_hitcall = if_else(assay_ac50 >= lower_bnd_um, 0, 1)) %>%
    select(casrn, assay_name, tested, cyto_hitcall)

  write_tsv(
    ac50_mat_long,
    here(Inter, "toxcast_activity_all_assays_invitrodb_v3-3.tsv")
  )

  write_tsv(
    lwr_bnd_hitc,
    here(Inter, "lower_cyto_pt_adjusted_hitcalls_all_assays_invitrodb_v3-3.tsv")
  )
}


# Identify all cytotoxicity assays ---------------------------     ---------------------------     ---------------------------
## Not just the 87 used to calculate the cytotoxicity point


## List the fields in the assay_component_endpoint table
tcplListFlds("assay_component_endpoint")

## Get aeids for assays whose intented_target_family_sub == cytotoxicity
## Should retrieve 117 assays
cyto_aeids <- tcplLoadAeid(
  fld = "intended_target_family_sub",
  val = "cytotoxicity"
)

## Get aeids for assays used to calculate burst (cytotoxicity point)
## Should retrieve 87 assays
burst_aeids <- tcplLoadAeid(fld = "burst_assay", val = 1)


# Retrieve information not available via tcpl ----------------------------------------------------------


## Connect to invitroDB via MariaDB
con <- DBI::dbConnect(MariaDB(),
                      dbname = dbconfig[['database']],
                      host = dbconfig[['host']],
                      port = dbconfig[['dbport']],
                      user = dbconfig[['user']],
                      password = dbconfig[['password']]
)

caution_assays_query <- '
SELECT
    assay_component_endpoint.assay_component_endpoint_name,
    assay_component_endpoint.assay_component_endpoint_desc
FROM assay_component_endpoint
WHERE assay_component_endpoint.assay_component_endpoint_desc LIKE "%Use data with caution.%"
'

caution_assays <- dbGetQuery(con, caution_assays_query)

dbDisconnect(con)
