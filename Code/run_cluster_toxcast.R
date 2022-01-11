# Will convert this to a Notebook once I have Pandoc working on the server

# Initialise functions in swe_functions.R ---------------------------     ---------------------------     ---------------------------

# Initializes functions and provides the following variables:
# catmos_fullset ac50_mat_long lwr_bnd_hitc
source("Code/swe_functions.R")

# Shouldn't be needed after the swe_function call, but just in case...
if(!(exists("ac50_mat_long") & exists("lwr_bnd_hitc"))){
  source("Code/retrieve_toxcast_data.R")
}

## Identifies ToxCast assays enriched for activity in (very) acutely toxic chemicals
if(!(exists("toxic_enriched_assays") & exists("very_toxic_enriched_assays"))){
  source(here("Code/enriched_assays_and_pathways.R"))
}

# Repeated from retrieve_toxcast_data.R for reference.
# This is called from swe_functions.R to load ToxCast data.
Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


# Load data files ---------------------------     ---------------------------     ---------------------------


## If cluster assignment file doesn't exist run catmos_chem_clustering.R
if (!file.exists(here(Output, "Toxprint/Cluster_assignment.tsv"))) {
  source(here("Code/catmos_chem_clustering.R"))
} else {
  chem_clusters <- read_tsv(
    here(Output, "Toxprint/Cluster_assignment.tsv"),
    col_types = cols()
  )
}

## If ToxCast assays enriched in each chemical cluster has been generated run cluster_enrichments.R
if (!file.exists(here(Output, "toxcast_assays_enriched_in_chemical_clusters_cut_height_0.7.tsv"))) {
  source(here("Code/cluster_enrichments.R"))
} else {
  cluster_tc_enrichments <- read_tsv(
    here(Output, "toxcast_assays_enriched_in_chemical_clusters_cut_height_0.7.tsv"),
    col_types = cols()
  )
}


# Process all assays ---------------------------     ---------------------------     ---------------------------


toxic_tc_fdr <- toxic_enriched_assays %>%
  full_join(very_toxic_enriched_assays, by = "chemotype") %>%
  mutate(tox_fdr = FP.x/(TP.x+FP.x), vt_fdr = FP.y/(TP.y+FP.y)) %>%
  select(chemotype, tox_fdr, vt_fdr)

filtered_cluster_tc_enrichments <- cluster_tc_enrichments %>%
  filter(mcc > 0.1 & CT_total > 3)

output_dir <- "OutputFiles/AllAssays"

source(here("Code/cluster_toxcast.R"))
length(intersect(minimal_assays, unique_mc_assays))
setdiff(minimal_assays, unique_mc_assays)

toxic_tc_fdr_count
toxcast_chems
toxcast_chems2
chem_tox_activity_summary
minimal_cluster_assays_summary
dtxsid_mismatches
chnm_mismatches
# Should match toxcast_chems2
cluster_activity_short_summary


write_tsv(cluster_activity_summary, here(Output, "cluster_activity_summary_all_assays.tsv"))


# Repeat excluding cytotox and caution assays ---------------------------     ---------------------------     ---------------------------


exclude_assays <- union(cyto_aeids$aenm, caution_assays$assay_component_endpoint_name)

toxic_tc_fdr <- toxic_enriched_assays %>%
  full_join(very_toxic_enriched_assays, by = "chemotype") %>%
  filter(!(chemotype %in% exclude_assays)) %>%
  mutate(tox_fdr = FP.x/(TP.x+FP.x), vt_fdr = FP.y/(TP.y+FP.y)) %>%
  select(chemotype, tox_fdr, vt_fdr)

filtered_cluster_tc_enrichments <- cluster_tc_enrichments %>%
  filter(mcc > 0.1 & CT_total > 3 & !(chemotype %in% exclude_assays))

lwr_bnd_hitc <- lwr_bnd_hitc %>%
  filter(!(assay_name %in% exclude_assays))

output_dir <- "OutputFiles/FilteredAssays"

source(here("Code/cluster_toxcast.R"))
length(intersect(minimal_assays, unique_mc_assays))
setdiff(minimal_assays, unique_mc_assays)

toxic_tc_fdr_count
toxcast_chems
toxcast_chems2
chem_tox_activity_summary
minimal_cluster_assays_summary
dtxsid_mismatches
chnm_mismatches
# Should match toxcast_chems2
cluster_activity_short_summary

## Write out summary files needed by "create_figures_2_and_3.R"
write_tsv(toxcast_chems, here(Output, "toxcast_chems_filtered_assays.tsv"))
write_tsv(chem_tox_activity_summary, here(Output, "chem_tox_activity_summary_filtered_assays.tsv"))
saveRDS(minimal_assays, here(Output, "minimal_assays_without_using_structural_clusters_filtered_assays.rds"))
saveRDS(unique_mc_assays, here(Output, "minimal_assays_with_using_structural_clusters_filtered_assays.rds"))
write_tsv(minimal_cluster_assays_summary, here(Output, "minimal_cluster_assays_summary_filtered_assays.tsv"))
write_tsv(cluster_activity_summary, here(Output, "cluster_activity_summary_filtered_assays.tsv"))

cluster_activity_summary %>%
  ungroup() %>%
  mutate(DSSTox_SID = if_else(is.na(dsstox_substance_id), DTXSID, dsstox_substance_id),
         chemical_name = if_else(is.na(chnm), PREFERRED_NAME, chnm),
         cluster_id = cluster_id_0.7) %>%
  select(DSSTox_SID, chemical_name, CASRN, very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category, toxic, cluster_id, assay) %>%
  write_tsv(here(Output, "chemical_cluster_assay_summary.tsv"))


# Examples ---------------------------     ---------------------------     ---------------------------


# Overall summaries
minimal_cluster_assays %>%
  group_by(assay) %>%
  summarise(n = n_distinct(cluster_id)) %>%
  arrange(desc(n))

minimal_cluster_assays %>%
  inner_join(chem_cluster_tox_activity) %>%
  group_by(cluster_id, assay) %>%
  summarise(cluster_assay_TP = sum(TP), cluster_assay_FP = sum(FP)) %>%
  group_by(cluster_id) %>%
  summarise(m = max(cluster_assay_FP)) %>%
  filter(m > 0) %>%
  arrange(desc(m))

minimal_cluster_assays %>%
  group_by(cluster_id) %>%
  count() %>%
  filter(n > 2) %>%
  arrange(desc(n))

# Cluster summaries
cluster <- 1964

# View hit calls
minimal_cluster_assays %>%
  filter(cluster_id == cluster) %>%
  select(assay)

dir.create("OutputFiles/Figures")
png(filename = "OutputFiles/Figures/figure11.png", width = 800, height = 800)
x <- chem_assay_heatmap(filter(chem_cluster_tox_activity, cluster_id == cluster), hm_title = glue("Enriched Assays - Cluster {cluster}"),
                        axis_margins = c(15,30), dendrogram = "none")
dev.off()

png(filename = "OutputFiles/Figures/figure5C.png", width = 800, height = 900)
y <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id == cluster), hm_title = glue("All Activity - Cluster {cluster}"),
                        hm_breaks = c(-2,-1,0,1,2), hm_colors = c("grey90", "white", "lightblue2", "blue"),
                        axis_margins = c(10,30), dendrogram = "none")
dev.off()

#Custom code for vitamin-D clusters
assays <- c("ATG_PXR_TRANS_up","ATG_PXRE_CIS_up","ATG_VDR_TRANS_up","CEETOX_H295R_OHPREG_dn",
            "CEETOX_H295R_OHPROG_dn","NHEERL_MED_hDIO3_dn","NVS_NR_mERa","TOX21_DT40",
            "TOX21_DT40_100","TOX21_DT40_657","TOX21_SSH_3T3_GLI3_Antagonist",
            "TOX21_VDR_BLA_agonist_ch2","TOX21_VDR_BLA_agonist_ratio","TOX21_VDR_BLA_Antagonist_ch1")
cluster <- 1798
png(filename = "OutputFiles/Figures/figure8A.png", width = 800, height = 600)
x <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id == cluster & assay %in% assays),
                        hm_title = glue("All Activity - Cluster {cluster}"),
                        hm_breaks = c(-2,-1,0,1,2), hm_colors = c("grey90", "white", "lightblue2", "blue"),
                        axis_margins = c(10,30), dendrogram = "none")
dev.off()

cluster <- 1964
png(filename = "OutputFiles/Figures/figure8B.png", width = 800, height = 600)
y <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id == cluster & assay %in% assays),
                        hm_title = glue("All Activity - Cluster {cluster}"),
                        hm_breaks = c(-2,-1,0,1,2), hm_colors = c("grey90", "white", "lightblue2", "blue"),
                        axis_margins = c(17,30), dendrogram = "none")
dev.off()

png(filename = "OutputFiles/Figures/figure8.png", width = 800, height = 400)
x <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id %in% c(1798,1964) & assay %in% assays),
                        hm_title = glue("All Activity - Clusters 1798, 1964"),
                        hm_breaks = c(-2,-1,0,1,2), hm_colors = c("grey90", "white", "lightblue2", "blue"),
                        axis_margins = c(10,30), dendrogram = "none")
dev.off()


minimal_cluster_assays %>%
  inner_join(chem_cluster_tox_activity, by = c("cluster_id", "assay")) %>%
  filter(cluster_id == cluster & !nontoxic & activity >= 0) %>%
  group_by(cluster_id, CASRN, PREFERRED_NAME, nontoxic, toxic, EPA_category, GHS_category, activity) %>%
  summarise(n_distinct(assay))


minimal_cluster_assays %>%
  inner_join(chem_cluster_tox_activity, by = c("cluster_id", "assay")) %>%
  filter(cluster_id == cluster)

chem_cluster_tox_activity %>%
  filter(cluster_id == cluster & !nontoxic)

chem_clusters %>%
  filter(cluster_id_0.7 == cluster & !nontoxic)

cluster_activity_summary %>%
  ungroup() %>%
  filter(toxcast_cyto_hitcall == 1 & is.na(assay))

