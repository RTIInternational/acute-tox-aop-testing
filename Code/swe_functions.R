if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  glue,
  gplots,
  RColorBrewer,
  reshape2,
  tidyverse
)
options(dplyr.summarise.inform = FALSE)

if(!(exists("acute_proc"))){
  source("Code/retrieve_toxcast_data.R")
}

## Load global datasets for functions
#catmos_fullset <- read_delim("InputData/full_catmos_with_chem_names.tsv", "\t")
catmos_fullset <- acute_proc %>%
  left_join(chem_info, by = c("CASRN" = "casn")) %>%
  select(dsstox_substance_id, chnm, CASRN, very_toxic, nontoxic, LD50_mgkg, EPA_category, GHS_category)

## Plotting functions

# Trim hit call matrix by number of hits per row then plot heatmap
# Returns the hit call matrix used for plotting
plot_chem_cluster_heatmap <- function(hcm, min_clusters = 2, axis_margins = c(10,25),
                                      hm_title = "Community Heatmap"){
  #hcm <- hit_call_matrix[,V(subgraph)$name]
  if ("matrix" %in% class(hcm)){
    small_mat <- hcm[apply(hcm, 1, sum, na.rm = TRUE) >= min_clusters,]
    if (dim(small_mat)[1] > 1){
      new_mat <- t(small_mat)
    }else{
      new_mat <- t(hcm)
    }
    heatmap.2(new_mat, breaks = c(-1,0,1), trace = 'none', col = c("grey90", "blue"),
              cexRow = 1.2, margins = axis_margins, main = hm_title)
  }else{
    new_mat <- t(as.matrix(hcm))
  }
  new_mat
  #Ordering the rows outside of the heatmap function to avoid the dendograms on the output
  #row.order <- hclust(dist(new_mat))$order
  #col.order <- hclust(dist(t(new_mat)))$order
  #heatmap(new_mat[row.order,col.order], Rowv = NA, Colv = NA, scale = "none",
  #        margins = axis_margins)
}

table2heatmap <- function(df, min_clusters = 2, axis_margins = c(10,25),
                          hm_title = "Community Heatmap"){
  if (dim(df)[2] == 3){
    #df <-
    hm_tbl <- df %>%
      #pivot_wider(names_from = chem_cluster, values_from = value, values_fill = 0)
      pivot_wider(names_from = colnames(df)[2], values_from = colnames(df)[3], values_fill = 0)
    hm_mat <- as.matrix(hm_tbl[,2:dim(hm_tbl)[2]], dimnames = hm_tbl[,1])
  }else{
    print("Expected 3 column table as follows: row names, column names, cell values")
  }
  small_mat <- hm_mat[apply(hm_mat, 1, sum, na.rm = TRUE) >= min_clusters,]
  if (dim(small_mat)[1] > 1){
    new_mat <- t(small_mat)
  }else{
    new_mat <- t(hm_mat)
  }
  heatmap.2(new_mat, trace = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25), #c("grey90", "blue"), #coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
            cexRow = 1.2, margins = axis_margins, main = hm_title)
  return(new_mat)
}
#x <- table2heatmap(df)

#df <- filter(chem_cluster_tox_activity, cluster_id == 411)
#df <- filter(chem_cluster_all_activity, cluster_id == 6) #452 53 1990
chem_assay_heatmap <- function(df, tox_criteria = "toxic", axis_margins = c(15,25), hm_title = "Cluster Heatmap",
                               hm_breaks = c(-2,-1,0,1), hm_colors = c("grey90", "white", "blue"), ...){
  if (all(c("CASRN", "assay", "activity", tox_criteria) %in% colnames(df))){
    active_assays <- df %>%
      group_by(assay) %>%
      summarise(assay_activity = suppressWarnings(max(activity, na.rm = TRUE))) %>%
      filter(assay_activity >= 1 & assay_activity <= 2)
    hm_tbl <- df %>%
      select(CASRN, all_of(tox_criteria), assay, activity) %>%
#      filter(assay %in% active_assays$assay) %>%
      pivot_wider(names_from = assay, values_from = activity, values_fill = -1)
      #pivot_wider(names_from = colnames(df)[2], values_from = colnames(df)[3], values_fill = 0)
    hm_mat <- t(as.matrix(hm_tbl[,active_assays$assay]))#, dimnames = list(row = paste("C", hm_tbl$CASRN, sep = ""))))
    colnames(hm_mat) <- hm_tbl$CASRN
  }else{
    stop(glue("Expected the following column names: CASRN, assay, activity, {tox_criteria}"))
  }

  num_assays <- dim(hm_mat)[1]
  num_chems <- dim(hm_mat)[2]
  if (num_assays > 0 & num_chems > 0){
    tox_chem <- factor(hm_tbl[[2]], levels = c("Yes", "No", "Unknown"))
    levels(tox_chem) <-  brewer.pal(5, "Pastel1") #length(levels(tox_chem))
    if (num_assays > 1 & num_chems > 1){
      if (num_assays > 50){print(num_assays)}
      heatmap.2(hm_mat, trace = 'none', breaks = hm_breaks, col = hm_colors, key = FALSE,
                #colorRampPalette(brewer.pal(8, "Blues"))(25), #c("grey90", "blue"), #coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
                ColSideColors = as.character(tox_chem), cexRow = 1.2, margins = axis_margins, main = hm_title, ...)
    }else{
      long_hm <- melt(hm_mat)
      if (num_assays > 1){
        chem_name <- unique(long_hm$Var2)
        if (length(chem_name) != 1){
          warning(glue("Error transforming heatmap matrix for bar chart: {hm_title}"))
        }else{
          print(ggplot(long_hm, aes(x = Var1, y = value)) +
              geom_col() +
              coord_flip() +
              labs(x = "ToxCast Assay", y = glue("Active for {chem_name[1]}"), title = hm_title)
          )
        }
      }else{
        assay_name <- unique(long_hm$Var1)
        if (length(assay_name) != 1){
          warning(glue("Error transforming heatmap matrix for bar chart: {hm_title}"))
        }else{
          print(ggplot(long_hm, aes(x = Var2, y = value)) +
            geom_col() +
            labs(x = "Chemical (CASRN)", y = glue("Active in {assay_name[1]}"), title = hm_title)
          )
        }
      }
    }
  }
  return(hm_mat)
}
#x <- chem_assay_heatmap(filter(chem_cluster_tox_activity, cluster_id == 6)) #28
#y <- chem_assay_heatmap(filter(chem_cluster_all_activity, cluster_id == 452))  #452 53 1990

## Data preparation

# Create a modified hit call matrix based on enrichment for toxic chemicals
calc_tox_enrichment <- function(hit_calls){
  # Pivot table if needed & Grab copy of input hit call matrix for QC below
  # Assumes 3 columns in long form ordered as so: CASRN, endpoints, hits
  if (dim(hit_calls)[2] == 3){
    #hit_calls <- hit_call_toxcast
    colnames(hit_calls) <- c("CASRN", "endpoint", "hit")
    phc_617 <- hit_calls
    hit_calls <- hit_calls %>%
      pivot_wider(names_from = endpoint, values_from = hit, values_fill = 0)
  }else{
    #hit_calls <- hit_call_617
    phc_617 <- pivot_longer(hit_calls, !CASRN, names_to = "endpoint", values_to = "hit")
    # Assume that NAs are not hits (no worse that an item not shopped because it was out of stock)
    hit_calls[is.na(hit_calls)] <- 0
  }

  ## Identify pathways, phenotypes, assays
  # Merge hit call matrix with information on toxicity
  ann_hc <- catmos_fullset %>%
    filter(!is.na(very_toxic) & !is.na(nontoxic)) %>%
    mutate(toxicity_cat = if_else(very_toxic==TRUE, 'very_toxic', if_else(nontoxic==TRUE, 'nontoxic', 'toxic'))) %>%
    inner_join(hit_calls, by = "CASRN")
  hc_ordered <- ann_hc[,10:dim(ann_hc)[2]]

  # Calcualte enrichment and extract key fields
  #fisher.test(annotated_hit_call_617$nontoxic, pull(annotated_hit_call_617[,10]))
  hc_sum <- apply(hc_ordered, 2, sum)
  hc_final <- hc_ordered[,hc_sum!=0]
  nontox_enrich <- apply(hc_final, 2, function(x) fisher.test(ann_hc$nontoxic, x))
  verytox_enrich <- apply(hc_final, 2, function(x) fisher.test(ann_hc$very_toxic, x))
  toxcat_enrich <- apply(hc_final, 2, function(x) fisher.test(ann_hc$toxicity_cat, x))
  nontox_pval <- sapply(nontox_enrich, function(x) x$p.value)
  nontox_OR <- sapply(nontox_enrich, function(x) x$estimate)
  verytox_pval <- sapply(verytox_enrich, function(x) x$p.value)
  verytox_OR <- sapply(verytox_enrich, function(x) x$estimate)
  toxcat_pval <- sapply(toxcat_enrich, function(x) x$p.value)

  # Create queryable summary of enrichment for plotting
  enrichment <- bind_rows(
    data.frame(category = 'nontox', pvals = nontox_pval, odds_ratio = nontox_OR),
    data.frame(category = 'verytox', pvals = verytox_pval, odds_ratio = verytox_OR),
    data.frame(category = 'toxcat', pvals = toxcat_pval, odds_ratio = NA)
    )

  # Compare the new hit call matrix to the old one
  phc_final <- bind_cols(ann_hc[,"CASRN"], hc_final) %>%
    pivot_longer(!CASRN, names_to = "endpoint", values_to = "hit")

  hc_compare <- phc_617 %>%
    full_join(phc_final, by = c("CASRN", "endpoint"))

  hc_compare %>%
    filter(hit.x != hit.y | (is.na(hit.y) & !(endpoint %in% names(hc_sum[hc_sum==0])))) %>%
    group_by(hit.x, hit.y) %>%
    summarise(n(), n_distinct(endpoint), n_distinct(CASRN)) %>%
    print()

  hc_compare_cas <- hc_compare %>%
    filter(hit.x != hit.y | (is.na(hit.y) & !(endpoint %in% names(hc_sum[hc_sum==0])))) %>%
    distinct(CASRN) %>%
    left_join(catmos_fullset, by = "CASRN")

  hc_compare_endpoint <- hc_compare %>%
    filter(hit.x == 1 & (hit.x != hit.y | is.na(hit.y)))

    return(list(hit_calls = hc_final, enrichment = enrichment, hc_rowi = ann_hc[,1:9],
              hc_coli = tibble(nontox_pval = nontox_pval, nontox_OR = nontox_OR,
                                verytox_pval = verytox_pval, verytox_OR = verytox_OR,
                                toxcat_pval),
              hc_compare_cas = hc_compare_cas, hc_compare_endpoint = hc_compare_endpoint
              )
         )
}
#enrichment_results <- calc_tox_enrichment(hit_call_617)

plot_stats <- function (community_rules){
  print(community_rules %>%
          group_by(RHS) %>%
          summarise(min_lift = min(lift), med_lift = median(lift), max_lift = max(lift)) %>%
          pivot_longer(cols = -RHS, names_to = "stat", values_to = "value") %>%
          ggplot(aes(value)) +
          geom_histogram() +
          facet_grid(stat ~ .) +
          labs(title = "Lift distribution for RHS")
  )
  print(community_rules %>%
          group_by(LHS) %>%
          summarise(min_lift = min(lift), med_lift = median(lift), max_lift = max(lift)) %>%
          pivot_longer(cols = -LHS, names_to = "stat", values_to = "value") %>%
          ggplot(aes(value)) +
          geom_histogram() +
          facet_wrap(~ stat) +
          labs(title = "Lift distribution for LHS")
  )

  print(community_rules %>%
          group_by(RHS) %>%
          summarise(count = n()) %>%
          ggplot(aes(count)) +
          geom_histogram(binwidth = 2) +
          labs(title = "Chem count distribution for RHS", x = "Chem count for rule")
  )
  print(community_rules %>%
          group_by(LHS) %>%
          summarise(count = n()) %>%
          ggplot(aes(count)) +
          geom_histogram(binwidth = 2) +
          labs(title = "Chem count distribution for LHS", x = "Chem count for rule")
  )
  return (community_rules)
}

create_netvis_matrices <- function(enrichment_results, row_filter = "nontoxic", row_inverse = FALSE,
                                   main_pval = 0.05, nontox_pval = 0.05, verytox_pval = 0.05,
                                   nontox_OR = 0.5, verytox_OR = 2){
  hcm <- enrichment_results[["hit_calls"]]
  hcr <- enrichment_results[["hc_rowi"]]
  hcc <- enrichment_results[["hc_coli"]]
  enrichment <- enrichment_results[["enrichment"]]

  i <- pull(hcr[,row_filter])
  if (row_inverse){
    i <- !i
  }
  print(table(i))
  j <- hcc$toxcat_pval < main_pval &
    ((hcc$nontox_pval < nontox_pval & hcc$nontox_OR < nontox_OR) |
       (hcc$verytox_pval < verytox_pval & hcc$verytox_OR > verytox_OR))
  #hc_filtered <- hc_final[ann_hc$nontoxic,i]
  #hc_filtered <- as.matrix(hc_final[!ann_hc$nontoxic,i])
  print(table(j))

  hc_filtered <- as.matrix(hcm[i,j])
  rownames(hc_filtered) <- hcr$CASRN[i]
  txn_matrix <- as.matrix(mutate_all(hcm[i,j], as.logical))
  rownames(txn_matrix) <- hcr$CASRN[i]
  return(list(hit_calls = hc_filtered, transactions = txn_matrix))
}
#netvis_matrices <- create_netvis_matrices(enrichment_results, row_inverse = TRUE)
#dim(netvis_matrices[["transactions"]])

split_rules <- function(rules_df){
  rules_df %>%
    separate(rules, c('LHS', 'RHS'), '} => {') %>%
    mutate(LHS = str_remove_all(LHS, '[\\{\\}]'),
           RHS = str_remove_all(RHS, '[\\{\\}]')) %>%
    return()
}

## Working with chemotypes
get_chemotypes <- function(tp_enrichment){
  max_OR <- max(tp_enrichment$OR[is.finite(tp_enrichment$OR)])*2

  tp_filtered <- tp_enrichment %>%
    filter(mcc > 0.1 & CT_total > 3 & p.value < 0.0001)# %>%
    #mutate(OR = if_else(is.infinite(OR), max_OR, OR)) %>%
    #mutate(std_OR = OR/max_OR)

  print(tp_filtered %>%
          ggplot(aes(log10(OR))) +
          geom_histogram(binwidth = 1)
  )

  print(tp_filtered %>%
          ggplot(aes(mcc)) +
          geom_histogram(binwidth = 3)
  )

  print(tp_filtered %>%
          ggplot(aes(TP)) +
          geom_histogram(binwidth = 3)
  )

  tp_final <-  tp_filtered %>%
    mutate(node_id = as.character(.[[1]])) %>%
    select(node_id, chemotype, mcc)

  return(tp_final)
}

## Working with assays
get_minimal_assays <- function(chem_assay){
  minimal_assays <- character()
  #minimal_assays <- data.frame()

  while(dim(chem_assay)[1] > 0 & sum(chem_assay$TP) > 0){
    top_assays <- chem_assay %>%
      group_by(assay, tox_fdr, vt_fdr) %>%
      summarise(sum_TP = sum(TP), sum_FP = sum(FP)) %>%
      filter(sum_TP > 0) %>%
      arrange(desc(sum_TP - sum_FP), desc(sum_TP), tox_fdr, vt_fdr)

    minimal_assays <- append(minimal_assays, top_assays[[1,1]])
    #minimal_assays <- rbind(minimal_assays, top_assays[1,])

    covered_chems <- chem_assay %>%
      filter(assay == top_assays[[1,1]])
    chem_assay <- chem_assay %>%
      filter(!(CASRN %in% covered_chems$CASRN))
  }
  return(minimal_assays)
}

