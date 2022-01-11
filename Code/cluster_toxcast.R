## Initial Summaries
## Initial Summaries
# Enrichment for toxicity only, no structure
toxic_tc_fdr_count <- toxic_tc_fdr %>%
  filter(!is.na(tox_fdr) & !is.na(vt_fdr)) %>%
  count()

toxic_tc_fdr %>%
  pivot_longer(c(tox_fdr, vt_fdr)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(name ~ .)

# Structure-based enrichment
cluster_tc_enrichment_summary <- filtered_cluster_tc_enrichments %>%
  group_by(cluster_id) %>%
  summarise(effect_count = n(), min_p = min(p.value), max_p = max(p.value))

cluster_tc_enrichment_summary %>%
  ggplot(aes(x = effect_count)) +
  geom_histogram(binwidth = 1)

enriched_clusters_chems <- chem_clusters %>%
  filter(cluster_id_0.7 %in% cluster_tc_enrichment_summary$cluster_id &
           !nontoxic) %>%
  summarise(cluster_count = n_distinct(cluster_id_0.7), chem_count = n_distinct(CASRN))

## Replicate Mark's summary and connect to structural clusters
toxcast_chems <- catmos_fullset %>%
  inner_join(ac50_mat_long, by = c('CASRN' = "casrn")) %>%
  filter(tested == 1) %>%
  group_by(nontoxic) %>%
  summarise(tested = n_distinct(CASRN)) %>%
  ungroup() %>%
  full_join(
    catmos_fullset %>%
      inner_join(ac50_mat_long, by = c('CASRN' = "casrn")) %>%
      filter(tested == 1 & assay_ac50 < 1000000) %>%
      group_by(nontoxic) %>%
      summarise(active = n_distinct(CASRN)) %>%
      ungroup(),
    by = "nontoxic") %>%
  full_join(
    catmos_fullset %>%
      inner_join(lwr_bnd_hitc, by = c('CASRN' = "casrn")) %>%
      filter(tested == 1 & cyto_hitcall == 1) %>%
      group_by(nontoxic) %>%
      summarise(cyto_hitcall = n_distinct(CASRN)) %>%
      ungroup(),
    by = "nontoxic")

# Create summary table and confirm that numbers match summaries above
tox_activity <- catmos_fullset %>%
  left_join(ac50_mat_long, by = c('CASRN' = "casrn")) %>%
  left_join(lwr_bnd_hitc, by = c('CASRN' = "casrn", "assay_name" = "assay_name", "tested" = "tested")) %>%
  mutate(toxic = if_else(is.na(nontoxic), "Unknown", if_else(nontoxic, "No", "Yes")),
         active = if_else(assay_ac50 < 1000000, 1, 0),
         tested = if_else(is.na(tested), 0, as.double(tested)),
         assay = assay_name, activity = if_else(tested == 1, active + cyto_hitcall, -1)
         )

tox_activity_summary <- tox_activity %>%
  group_by(dsstox_substance_id, chnm, CASRN,
           very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category, toxic) %>%
  summarise(tested_toxcast = max(tested), active_toxcast = max(active),
            toxcast_cyto_hitcall = max(if_else(is.na(cyto_hitcall), 0, cyto_hitcall)))

toxcast_chems2 <- tox_activity_summary %>%
  group_by(toxic, tested_toxcast, active_toxcast, toxcast_cyto_hitcall) %>%
  summarise(n_distinct(CASRN))

## Finding minimal assay sets
# Assemble all data into one data frame
chem_tox_activity <- catmos_fullset %>%
  inner_join(lwr_bnd_hitc, by = c("CASRN" = "casrn")) %>%
  left_join(toxic_tc_fdr, by = c("assay_name" = "chemotype")) %>%
  mutate(toxic = if_else(is.na(nontoxic), "Unknown", if_else(nontoxic, "No", "Yes")),
         assay = assay_name, activity = if_else(tested == 1, cyto_hitcall, -1),
         TP = if_else(!nontoxic & activity == 1, 1, 0), FP = if_else(nontoxic & activity == 1, 1, 0)) %>%
  select(CASRN, chnm, nontoxic, toxic, EPA_category, GHS_category,
         assay, activity, tox_fdr, vt_fdr, TP, FP)

# Should match values from toxcast_chems2
chem_tox_activity_summary <- chem_tox_activity %>%
  filter(activity > 0) %>%
  group_by(toxic) %>%
  summarise(n_distinct(CASRN)) %>%
  ungroup()

# Calculate the minimal number of ToxCast assays to capture all possible acute tox chemicals
minimal_assays <- chem_tox_activity %>%
  filter(activity == 1 & !is.na(TP)) %>%
  select(CASRN, assay, tox_fdr, vt_fdr, TP, FP) %>%
  get_minimal_assays()

# True positives vs. false positives for selected assays
# TP should match match values from toxcast_chems2 & chem_tox_activity_summary

chem_tox_activity_summary[chem_tox_activity_summary$toxic=="Yes","all_assay_positives"] <- chem_tox_activity %>%
  filter(TP == 1 & assay %in% minimal_assays) %>%
  summarise(n_distinct(CASRN))

chem_tox_activity_summary[chem_tox_activity_summary$toxic=="No","all_assay_positives"] <- chem_tox_activity %>%
  filter(FP == 1 & assay %in% minimal_assays) %>%
  summarise(n_distinct(CASRN))

#Repeat process for assays enriched for structural clusters
chem_cluster_all_activity <- tox_activity %>%
  inner_join(chem_clusters, by = c("CASRN", "very_toxic", "nontoxic", "LD50_mgkg", "EPA_category", "GHS_category")) %>%
  mutate(cluster_id = cluster_id_0.7) %>%
  select(cluster_id, CASRN, PREFERRED_NAME, nontoxic, toxic, EPA_category, GHS_category,
         assay, activity)

chem_cluster_tox_activity <- filtered_cluster_tc_enrichments %>%
  inner_join(chem_clusters, by = c("cluster_id" = "cluster_id_0.7")) %>%
  left_join(toxic_tc_fdr, by = "chemotype") %>%
  inner_join(lwr_bnd_hitc, by = c("CASRN" = "casrn", "chemotype" = "assay_name")) %>%
  mutate(toxic = if_else(is.na(nontoxic), "Unknown", if_else(nontoxic, "No", "Yes")),
         assay = chemotype, activity = if_else(tested == 1, cyto_hitcall, -1),
         TP = if_else(!nontoxic & activity == 1, 1, 0), FP = if_else(nontoxic & activity == 1, 1, 0)) %>%
  select(cluster_id, CASRN, PREFERRED_NAME, nontoxic, toxic, EPA_category, GHS_category,
         assay, activity, tox_fdr, vt_fdr, TP, FP)

chem_cluster_tox_activity %>%
  filter(activity > 0) %>%
  group_by(toxic) %>%
  summarise(n_distinct(CASRN))

# For testing the get_minimal_assays function
# cluster <- 861
# chem_assay <- chem_cluster_tox_activity %>%
#   mutate(TP = if_else(!nontoxic & activity == 1, 1, 0), FP = if_else(nontoxic & activity == 1, 1, 0)) %>%
#   filter(cluster_id == cluster & activity == 1) %>%
#   select(cluster_id, CASRN, assay, tox_fdr, vt_fdr, TP, FP)

# Find minimal assay set for each structural cluster
if (!dir.exists(glue("{output_dir}/chem_enriched_assay_heatmaps"))){
  dir.create(glue("{output_dir}/chem_enriched_assay_heatmaps"), recursive = TRUE, mode = "0750")
}
pdf(glue("{output_dir}/chem_enriched_assay_heatmaps.pdf"), 10, 8)
minimal_cluster_assay_list <- list()
for(cluster in sort(unique(chem_cluster_tox_activity$cluster_id))){
  x <- chem_assay_heatmap(filter(chem_cluster_tox_activity, cluster_id == cluster), hm_title = glue("Enriched Assay Heatmap for Cluster {cluster}"))
  write_tsv(as.data.frame(cbind(rownames(x), x)), here(glue("{output_dir}/chem_enriched_assay_heatmaps/cluster{cluster}.tsv")))
  ma <- tibble(cluster, chem_cluster_tox_activity %>%
     filter(cluster_id == cluster & activity == 1) %>%
     select(cluster_id, CASRN, assay, tox_fdr, vt_fdr, TP, FP) %>%
     get_minimal_assays()
  )

  minimal_cluster_assay_list[[cluster]] <- ma
}
minimal_cluster_assays <- bind_rows(minimal_cluster_assay_list)
colnames(minimal_cluster_assays) <- c("cluster_id", "assay")
dev.off()

# Create heatmaps for all ToxCast activity for comparison with above
if (!dir.exists(glue("{output_dir}/chem_all_activity_heatmaps"))){
  dir.create(glue("{output_dir}/chem_all_activity_heatmaps"), recursive = TRUE, mode = "0750")
}
pdf(glue("{output_dir}/chem_all_activity_heatmaps.pdf"), 10, 8)
for(cluster in sort(unique(chem_cluster_all_activity$cluster_id))){
  y <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id == cluster), hm_title = glue("All Activity Heatmap for Cluster {cluster}"),
                          hm_breaks = c(-2,-1,0,1,2), hm_colors = c("grey90", "white", "lightblue2", "blue"))
  write_tsv(as.data.frame(cbind(rownames(y), y)), here(glue("{output_dir}/chem_all_activity_heatmaps/cluster{cluster}.tsv")))
}
dev.off()

# Summarize minimal assays by cluster
unique_mc_assays <- unique(minimal_cluster_assays$assay)

minimal_cluster_assays %>%
  group_by(cluster_id) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram() +
  labs(x = "Number of assays per cluster", y = "Number of clusters")

minimal_cluster_assays_summary <- minimal_cluster_assays %>%
  group_by(cluster_id) %>%
  count() %>%
  group_by(n) %>%
  count()

chem_tox_activity_summary[chem_tox_activity_summary$toxic=="Yes","cluster_assay_positives"] <- minimal_cluster_assays %>%
  inner_join(chem_cluster_tox_activity, by = c("cluster_id", "assay")) %>%
  filter(TP == 1) %>%
  summarise(n_distinct(CASRN))

chem_tox_activity_summary[chem_tox_activity_summary$toxic=="No","cluster_assay_positives"] <- minimal_cluster_assays %>%
  inner_join(chem_cluster_tox_activity, by = c("cluster_id", "assay")) %>%
  filter(FP == 1) %>%
  summarise(n_distinct(CASRN))

# Expand summary table to include clusters and minimal assays
cluster_activity_summary <- tox_activity_summary %>%
  left_join(chem_clusters, by = c("CASRN", "very_toxic", "nontoxic", "LD50_mgkg", "EPA_category", "GHS_category")) %>%
  left_join(minimal_cluster_assays, by = c("cluster_id_0.7" = "cluster_id")) %>%
  select(dsstox_substance_id, DTXSID, chnm, PREFERRED_NAME, CASRN, very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category,
         toxic, tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7, assay)

# cluster_activity_summary <- catmos_fullset %>%
#   inner_join(chem_clusters %>% select(CASRN, cluster_id_0.7), by = c("CASRN")) %>%
#   inner_join(filtered_cluster_tc_enrichments, by = c("cluster_id_0.7" = "cluster_id")) %>%
#   left_join(ac50_mat_long, by = c("CASRN" = "casrn", "chemotype" = "assay_name")) %>%
#   inner_join(lwr_bnd_hitc, by = c("CASRN" = "casrn", "chemotype" = "assay_name", "tested")) %>%
#   mutate(
#     active_toxcast = if_else(assay_ac50 < 1000000, 1, 0),
#     tested = if_else(is.na(tested), 0, as.double(tested))
#   ) %>%
#   rename(
#     assay = chemotype,
#     tested_toxcast = tested,
#     toxcast_cyto_hitcall = cyto_hitcall
#   ) %>%
#   select(
#     cluster_id_0.7, dsstox_substance_id, DTXSID, chnm, PREFERRED_NAME, CASRN,
#     very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category,
#     assay, tested_toxcast, active_toxcast, toxcast_cyto_hitcall
#   )

dtxsid_mismatches <- filter(cluster_activity_summary, dsstox_substance_id != DTXSID)
chnm_mismatches <- filter(cluster_activity_summary, chnm != PREFERRED_NAME)
length(unique(cluster_activity_summary$CASRN))
summarize(catmos_fullset, n_distinct(CASRN))

# Should match toxcast_chems2
cluster_activity_short_summary <- cluster_activity_summary %>%
  group_by(toxic, tested_toxcast, active_toxcast, toxcast_cyto_hitcall) %>%
  summarise(n_distinct(CASRN))
