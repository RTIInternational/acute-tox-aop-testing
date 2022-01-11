# TODO Save relevant datasets in cluster_toxcast.R and read them to save re-processing
# variables used: minimal_cluster_assays cluster_activity_summary
# Add code to run retrieve_toxcast_data.R for the following variables:
# cyto_aeids  ac50_mat_long

# Cytotoxicity assays in minimal sets
# Should return no results when run on the filtered assay results
minimal_cluster_assays %>%
  inner_join(cyto_aeids, by = c("assay" = "aenm")) %>%
  summarise(n_distinct(cluster_id), n_distinct(assay))

minimal_cluster_assays %>%
  left_join(cyto_aeids, by = c("assay" = "aenm")) %>%
  mutate(cytotox_assay = if_else(is.na(intended_target_family_sub), 0, 1)) %>%
  group_by(cluster_id) %>%
  summarise(assay_count = n(), cytotox_assay_count = sum(cytotox_assay)) %>%
  filter(assay_count == cytotox_assay_count)

# Calculate cytotox score for chemicals
cytotox_activity <- cluster_activity_summary %>%
  inner_join(ac50_mat_long, by = c("CASRN" = "casrn")) %>%
  filter(assay_name %in% cyto_aeids$aenm & tested == 1) %>%
  group_by(dsstox_substance_id, DTXSID, chnm, PREFERRED_NAME, CASRN,
           very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category,
           tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7) %>%
  summarise(cyto_ac50 = median(assay_ac50))

# Summarize activity in minimal assays
x <- cluster_activity_summary %>%
  inner_join(ac50_mat_long, by = c("CASRN" = "casrn", "assay" = "assay_name")) %>%
  ungroup() %>%
  filter(tested == 1 & assay_ac50 < 1000000)

x %>%
  group_by(nontoxic) %>%
  summarize(n_distinct(CASRN))

x %>%
  group_by(GHS_category) %>%
  summarize(n_distinct(CASRN), median(assay_ac50))

x %>%
  ggplot(mapping = aes(x = as.factor(GHS_category), y = assay_ac50)) +
  geom_boxplot(notch = TRUE)

x2 <- cluster_activity_summary %>%
  inner_join(ac50_mat_long, by = c("CASRN" = "casrn", "assay" = "assay_name")) %>%
  ungroup() %>%
  filter(tested == 1 & assay_ac50 >= 1000000)

x2 %>%
  group_by(nontoxic) %>%
  summarize(n_distinct(CASRN))

x2 %>%
  group_by(GHS_category) %>%
  summarize(n_distinct(CASRN), median(assay_ac50))

# x %>%
#   ggplot(mapping = aes(x = GHS_category, y = assay_ac50)) +
#   geom_point()

# Summarize activity in cytotox assays
y <- cytotox_activity %>%
  ungroup() %>%
  filter(cyto_ac50 < 1000000)

y %>%
  group_by(nontoxic) %>%
  summarize(n_distinct(CASRN))

y %>%
  group_by(GHS_category) %>%
  summarize(n_distinct(CASRN), median(cyto_ac50))

y %>%
  ggplot(mapping = aes(x = as.factor(GHS_category), y = cyto_ac50)) +
  geom_boxplot()

x %>%
  ggplot(mapping = aes(x = GHS_category, y = assay_ac50)) +
  geom_point()

# Summarize activity
# old_fas <- final_activity_summary
final_activity_summary <- cluster_activity_summary %>%
  ungroup() %>%
  left_join(filter(ac50_mat_long, tested == 1), by = c("CASRN" = "casrn", "assay" = "assay_name")) %>%
  group_by(dsstox_substance_id, DTXSID, chnm, PREFERRED_NAME, CASRN,
           very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category, toxic,
           tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7) %>%
  summarise(noncyto_tox_ld50 = suppressWarnings(min(assay_ac50, na.rm = TRUE))) %>%
  left_join(cyto_pt, by = c("CASRN" = "casn")) %>%
  mutate(final_ac50 = if_else(is.na(cyto_pt_um) & noncyto_tox_ld50 < 1000, noncyto_tox_ld50,
                                if_else(noncyto_tox_ld50 < cyto_pt_um, noncyto_tox_ld50, cyto_pt_um)),
         final_activity = if_else(final_ac50 < 1000, 1, 0)) %>%
  ungroup()
#anti_join(old_fas, final_activity_summary) %>% select(CASRN, tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7, noncyto_tox_ld50, cyto_pt_um, lower_bnd_um, final_ac50, final_activity)
#anti_join(final_activity_summary, old_fas) %>% select(CASRN, tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7, noncyto_tox_ld50, cyto_pt_um, lower_bnd_um, final_ac50, final_activity)

# Chemicals have activity in ToxCast but not in assays selected for their cluster
# These also have no cytotox value, so the final activity is NA (See query below)
check_missing_actives <- final_activity_summary %>%
  filter(is.na(final_ac50) & active_toxcast == 1)

final_activity_summary %>%
  filter(cluster_id_0.7 == 562) %>% # 1407, 506, 562, 246
  select(CASRN, tested_toxcast, active_toxcast, toxcast_cyto_hitcall, cluster_id_0.7, noncyto_tox_ld50, cyto_pt_um, lower_bnd_um, final_ac50, final_activity)

# Confirming that the counts from the final table match those from the previous steps
final_activity_short_summary <- final_activity_summary %>%
  group_by(toxic, tested_toxcast, active_toxcast, toxcast_cyto_hitcall) %>%
  summarise(n(), n_distinct(CASRN))

# 7997 chemicals not tested + 38 chemicals not tested in appropriate assay
# 3957 chemicals evaluated 2028 with activity in tox assays, 1929 with no activity
final_activity_summary %>%
  group_by(tested_toxcast, final_activity) %>%
  summarize(n(), n_distinct(CASRN))

# Comparison of activity with toxicity
final_activity_summary %>%
  filter(!is.na(final_activity)) %>%
  group_by(toxic, final_activity) %>%
  summarize(n(), n_distinct(CASRN))

x <- matrix(c(1406, 649, 621, 1279), nrow = 2,
            dimnames = list(ToxCast = c("Active", "NotActive"),
                            ATWG_class = c("Tox", "NonTox")))
fisher.test(x, alternative = "greater")
fisher.test(x, simulate.p.value=TRUE, B = 20000, alternative = "greater")

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  group_by(toxic, final_activity) %>%
  summarize(num_chems = n()) %>%
  pivot_wider(names_from = toxic, values_from = num_chems) %>%
  select(-final_activity) %>%
  as.matrix() %>%
  fisher.test(alternative = "greater")


final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  group_by(GHS_category, final_activity) %>%
  summarize(num_chems = n()) %>%
  pivot_wider(names_from = GHS_category, values_from = num_chems) %>%
  select(-final_activity) %>%
  as.matrix() %>%
  fisher.test(simulate.p.value=TRUE, B = 20000, alternative = "greater")

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  group_by(GHS_category) %>%
  summarize(n_distinct(CASRN), median(final_ac50))

fas_plot <- final_activity_summary %>%
  filter(!is.na(final_activity) & !is.na(GHS_category) & toxic != "Unknown") %>%
  ggplot(mapping = aes(x = as.factor(GHS_category), y = final_ac50)) +
  geom_boxplot()

png(filename = "OutputFiles/Figures/figure12A.png", width = 480, height = 280)
fas_plot +
  scale_y_continuous(trans='log10') +
  labs(x = "GHS Category", y = "log10 AC50")
dev.off()

fas_plot +
  geom_point(aes(x = as.factor(GHS_category), y = median(final_ac50)))

final_activity_summary %>%
  filter(!is.na(final_activity) & !is.na(GHS_category) & toxic != "Unknown") %>%
  ggplot(mapping = aes(x = as.factor(GHS_category), y = final_ac50)) +
  geom_violin(adjust = 4)

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  ggplot(mapping = aes(x = as.factor(GHS_category), y = final_ac50)) +
  geom_dotplot(method = "histodot", binwidth = 10, binaxis = "y", width = 50, stackdir = "center", dotsize = 1/5)

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  ggplot(mapping = aes(x = final_ac50, fill = as.factor(GHS_category))) +
  geom_histogram(binwidth = 5)

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown" & final_ac50 < 100) %>%
  ggplot(mapping = aes(x = final_ac50, fill = as.factor(GHS_category))) +
  geom_histogram(binwidth = 2)

final_activity_summary %>%
  filter(!is.na(final_activity) & toxic != "Unknown") %>%
  ggplot(mapping = aes(x = final_ac50, fill = as.factor(GHS_category))) +
  geom_bar(position = "fill") +
  scale_x_binned(n.breaks = 100)

# Compare structural category with tox categories
ghs_clusters <- final_activity_summary %>%
  filter(!is.na(final_activity) & !is.na(GHS_category) & toxic != "Unknown") %>%
  group_by(cluster_id_0.7, GHS_category) %>%
  summarize(num_chems = n()) %>%
  pivot_wider(names_from = GHS_category, values_from = num_chems, values_fill = 0) %>%
  ungroup()

ghs_cluster_mat <- as.matrix(select(ghs_clusters, c("1", "2", "3", "4", "5")))
rownames(ghs_cluster_mat) <- paste("C", ghs_clusters$cluster_id_0.7, sep = "-")
fisher.test(ghs_cluster_mat, simulate.p.value=TRUE, B = 20000, alternative = "greater")

# Used to create figure12B.png (file exported from RStudio)
cluster_total <- apply(ghs_cluster_mat, 1, sum)
ghs_cluster_percentage <- ghs_cluster_mat/cluster_total
row_index <- order(ghs_cluster_percentage[,1], ghs_cluster_percentage[,2], ghs_cluster_percentage[,3], ghs_cluster_percentage[,4], ghs_cluster_percentage[,5])
heatmap.2(ghs_cluster_percentage[row_index,],
          trace = 'none', Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          margins = c(2, 1), keysize = 1, key.title = NA)

# Summarize coverage of chemicals by structural clusters
tox_clusters <- final_activity_summary %>%
  filter(toxic == "Yes" & !is.na(cluster_id_0.7)) %>%
  distinct(cluster_id_0.7) %>%
  pull(cluster_id_0.7)

active_clusters <- final_activity_summary %>%
  filter(final_activity == 1 & !is.na(cluster_id_0.7)) %>%
  distinct(cluster_id_0.7) %>%
  pull(cluster_id_0.7)

active_tox_clusters <- intersect(tox_clusters, active_clusters)

final_activity_summary %>%
  filter(toxic == "Yes") %>%
  count()

final_activity_summary %>%
  filter(toxic == "Yes" & cluster_id_0.7 %in% tox_clusters) %>%
  count()
6299/6845
length(tox_clusters)

final_activity_summary %>%
  filter(toxic == "Yes" & cluster_id_0.7 %in% active_tox_clusters) %>%
  count()
3723/6845
length(active_tox_clusters)

#Check distribution of GHS categories within clusters
cluster_ghs_summary <- final_activity_summary %>%
  filter(!is.na(GHS_category) & cluster_id_0.7 %in% tox_clusters) %>% #toxic == "Yes" &
  group_by(cluster_id_0.7) %>%
  summarize(num_chems = n(), num_ghs_cats = n_distinct(GHS_category),
            min_ghs_cat = min(GHS_category), max_ghs_cat = max(GHS_category)) %>%
  ungroup()

cluster_ghs_summary %>%
  ggplot(aes(x = num_ghs_cats)) +
  geom_histogram(binwidth = 1)

cluster_ghs_summary %>%
  group_by(num_ghs_cats) %>%
  count()
1039/sum(1039,579,169,18) #toxic == "Yes"
(1039+579)/sum(1039,579,169,18) #toxic == "Yes"
274/sum(274,1046,405,74,10) #all chems in tox clusters
sum(274,1046)/sum(274,1046,405,74,10) #all chems in tox clusters

cluster_ghs_summary %>%
  mutate(ghs_spread = max_ghs_cat - min_ghs_cat) %>%
  group_by(ghs_spread) %>%
  count()
1039/sum(1039,489,224,53) #toxic == "Yes"
(1039+489)/sum(1039,489,224,53) #toxic == "Yes"

cluster_ghs_summary %>%
  filter(min_ghs_cat > 3) %>%
  mutate(ghs_spread = max_ghs_cat - min_ghs_cat) %>%
  group_by(ghs_spread) %>%
  count()

x <- cluster_ghs_summary %>%
  filter(min_ghs_cat > 3)

cluster_ghs_summary %>%
#  filter(num_ghs_cats == 1 & max_ghs_cat - min_ghs_cat > 0)
  filter(is.na(max_ghs_cat)) %>%
  group_by(num_ghs_cats) %>%
  count()
