# ---------------------- Hypergeometric Enrichment Test ----------------------
hypergeometric_test <- function(df, cell_type, brain_region) {
  N <- sum(df)  # total count
  K <- sum(df[cell_type, ])
  n <- sum(df[, brain_region])
  k <- df[cell_type, brain_region]
  if (k == 0) return(1)
  phyper(k - 1, K, N - K, n, lower.tail = FALSE)
}

# Compute enrichment p-values across all cell type and region combinations
results <- expand.grid(Cell_Type = rownames(df_raw), Brain_Region = colnames(df_raw)) %>%
  rowwise() %>%
  mutate(P_Value = hypergeometric_test(df_raw, Cell_Type, Brain_Region)) %>%
  ungroup() %>%
  mutate(
    Adj.P.value = p.adjust(P_Value, method = "BH"),
    log.adj.p = -log10(pmax(Adj.P.value, 1e-300))
  )

# Significant enriched regions
sig_res <- filter(results, Adj.P.value < 0.05)
roi_regions <- unique(sig_res$Brain_Region)

# ---------------------- Permutation Test for Enrichment ----------------------
permutation_test <- function(region, data, n_perm = 1000) {
  region_data <- filter(data, brain_region == region)
  all_slices <- unique(region_data$slice)
  region_data <- filter(data, slice %in% all_slices)
  
  # Observed counts
  obs <- region_data %>%
    mutate(region_celltype = paste(brain_region, cell_type, sep = "_")) %>%
    count(region_celltype, name = "count")
  
  # Permutation background
  bg_list <- list()
  for (i in seq_len(n_perm)) {
    permuted <- region_data %>%
      mutate(
        cell_type = sample(cell_type),
        region_celltype = paste(brain_region, cell_type, sep = "_")
      ) %>%
      count(region_celltype, name = "count")
    
    for (j in seq_len(nrow(permuted))) {
      key <- permuted$region_celltype[j]
      bg_list[[key]] <- c(bg_list[[key]], permuted$count[j])
    }
  }
  
  # Compute empirical p-values
  obs$p_value <- sapply(obs$region_celltype, function(key) {
    bg_values <- bg_list[[key]]
    obs_val <- obs$count[obs$region_celltype == key]
    if (length(bg_values) > 0) {
      (sum(obs_val < bg_values) + 1) / (length(bg_values) + 1)  # smoothed p
    } else {
      1
    }
  })
  
  # Extract results for this region
  obs %>%
    filter(grepl(paste0("^", region, "_"), region_celltype)) %>%
    mutate(
      region = sub("_.*", "", region_celltype),
      cell_type = sub("^[^_]*_", "", region_celltype),
      log.p = -log10(pmax(p_value, 1e-300))
    ) %>%
    select(region_celltype, cell_type, region, p_value, log.p)
}

# Run permutation tests for significant regions
res_list <- lapply(roi_regions, permutation_test, data = pers)
result <- bind_rows(res_list) %>% arrange(cell_type, p_value)