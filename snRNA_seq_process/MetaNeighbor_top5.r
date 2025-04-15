suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
# suppressMessages(library(MetaNeighbor))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(Matrix))
suppressMessages(library(viridis))
suppressMessages(library(cowplot))
suppressMessages(library(ggsci))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
#suppressMessages(library(LSD))
suppressMessages(library(readxl))
suppressMessages(library(ggrepel))
suppressMessages(library(harmony))
suppressMessages(library(scrattch.hicat))
# suppressMessages(library(getopt))

reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
  
  conf1 <- table(cl, ref.cl) # a count table with cl as row and ref.cl as column
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/") # make the sum of each row = 1, the ratio of one cl related to each ref.cl
  conf2 <- reorder_matrix(conf1) # order the matrix to put high value to the diagonal
  
  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x) # 一个完全展开的2列data.frame
    min.prop <- apply(grid1, 1, min)
  })
  
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))
  # cl.prop.cocl.m 记录了cl中任意两类之间的共聚类的比例,取的是最小值之和,个人觉得这个计算方法很奇怪
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}


#' Identify a highly variable gene set
#'
#' Identifies genes with high variance compared to their median expression
#' (top quartile) within each experimentCertain function
#'
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix data
#' @param exp_labels character vector that denotes the source (Study ID) of
#' each sample.
#' @param min_recurrence Number of studies across which a gene must be detected
#' as highly variable to be kept. By default, only genes that are variable
#' across all studies are kept (intersection).
#' @param downsampling_size Downsample each study to downsampling_size
#' samples without replacement. If set to 0 or value exceeds dataset size,
#' no downsampling is applied.
#'
#' @return The output is a vector of gene names that are highly variable in
#' every experiment (intersect)
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' var_genes
#'
#' @export
#'

variableGenes <- function(dat, i = 1, exp_labels,
                          min_recurrence = length(unique(exp_labels)),
                          downsampling_size = 10000) {

    dat <- SummarizedExperiment::assay(dat, i = i)
    experiments <- unique(exp_labels)

    #check length of exp_labels equal # of samples
    if(length(exp_labels) != length(colnames(dat))){
        stop('experiment_labels length does not match number of samples')
    }
    if (min_recurrence > length(experiments)) {
        stop('min_recurrence should be smaller or equal to the number of datasets')
    }

    gene_stats <- list()
    for(exp in experiments){
        keep <- which(exp_labels == exp)
        if (downsampling_size > 0 & downsampling_size < length(keep)) {
            keep <- sample(keep, downsampling_size, replace = FALSE)
        }
        gene_stats[[exp]] <- variable_genes_single_exp(dat[, keep])
    }
    gene_stats <- dplyr::bind_rows(gene_stats, .id = "study_id")
    `%>%` <- dplyr::`%>%`
    gene_stats <- gene_stats %>%
        dplyr::group_by(gene) %>%
        dplyr::summarize(recurrence = sum(is_hvg),
                         score = mean(var_quant)) %>%
        dplyr::filter(recurrence >= min_recurrence) %>%
        dplyr::arrange(desc(recurrence), desc(score))
    result <- dplyr::pull(gene_stats, gene)
    return(result)
}

variable_genes_single_exp = function(data_subset) {
    data_subset <- as.matrix(data_subset)
    variance_data <- MN_rowVars(data_subset)
    median_data <- matrixStats::rowMedians(data_subset)
    quant_med <- unique(stats::quantile(
        median_data, probs = seq(0, 1, length = 11), type = 5
    ))
    
    # assign genes to 5 quantile bins based on median expression
    # remove bin with high expressing genes
    `%>%` <- dplyr::`%>%`
    result <- data.frame(gene = names(variance_data),
                         variance = variance_data,
                         bin_med = cut(median_data, c(-1,quant_med))) %>%
        dplyr::filter(bin_med != levels(bin_med)[length(levels(bin_med))]) %>%
        dplyr::group_by(bin_med) %>%
        dplyr::mutate(var_quant = (rank(variance)-1) / (length(variance)-1)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-variance, -bin_med) %>%
        dplyr::mutate(is_hvg = var_quant > 0.75)
    return(result)
}

MN_rowVars <- function(M) {
    if (methods::is(M, "dgCMatrix")) {
        M <- Matrix::t(M)
        result <- Matrix::colMeans(M**2) - Matrix::colMeans(M)**2
        result <- result * nrow(M) / (nrow(M)-1)
    } else {
        result <- matrixStats::rowVars(as.matrix(M))
        names(result) <- rownames(M)
    }
    return(result)
}


###loading function

makeClusterName <- function(study_id, cell_type) {
  if (length(study_id) != length(cell_type)) {
      stop("study_id and cell_type must have identical length!")
  }
  return(paste(standardizeLabel(study_id),
               standardizeLabel(cell_type),
               sep = "|"))
}

standardizeLabel <- function(labels, replace = "|", with = ".") {
  return(gsub(replace, with, labels, fixed = TRUE))
}

getStudyId <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 1))
}

getCellType <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "|", fixed = TRUE), "[", 2))
}

# Scale matrix such that all colums sum to 0 and have l2-norm of 1
normalize_cols <- function(M, ranked = TRUE) {
  result <- as.matrix(M)
  if (ranked) {
    result <- matrixStats::colRanks(result, ties.method = "average",
                                    preserveShape = TRUE)
  }
  result <- scale_cols(result)
  dimnames(result) <- dimnames(M)
  return(result)
}

scale_cols <- function(M) {
    cm <- colMeans(M)
    cnorm <- 1 / sqrt(colSums(M**2) - nrow(M) * cm**2)
    matrixStats::t_tx_OP_y(matrixStats::t_tx_OP_y(M, cm, "-"), cnorm, "*")
}

compute_aurocs <- function(votes, candidate_id = NULL) {
  if (is.null(candidate_id)) {
    positives <- design_matrix(rownames(votes))
  } else {
    positives <- as.matrix(candidate_id)
  }
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- crossprod(
    positives,
    matrixStats::colRanks(votes, ties.method = "average", preserveShape = TRUE)
  )
  colnames(sum_of_positive_ranks) <- colnames(votes)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

# Transform a vector with cell_type labels into a binary matrix
design_matrix <- function(cell_type) {
  cell_type <- as.factor(cell_type)
  if (length(levels(cell_type)) > 1) {
    result <- stats::model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- levels(cell_type)
  return(result)
}

# for (i in seq_len(n_candidates)) {
# for (j in seq(2, n_candidates)) {
#     best <- candidates[i]
#     contender <- candidates[j]
#     votes_best <- votes[names(votes) == best]
#     votes_contender <- votes[names(votes) == contender]
# #     if (length(votes_best) == 0 | length(votes_contender) == 0) {
# #         next
# #       }
#     auroc <- compute_aurocs(
#         as.matrix(c(votes_best, votes_contender)),
#         as.matrix(rep(c(1, 0), c(length(votes_best), length(votes_contender))))
#     )
#     auroc_matrix[best, contender] <- auroc
#     auroc_matrix[contender, best] <- 1 - auroc  # 对称矩阵
# }
# }

###

compute_1v1_aurocs <- function(votes, aurocs) {
  result <- matrix(NA, nrow(aurocs), ncol(aurocs), dimnames = dimnames(aurocs))
  for (i in seq_len(ncol(aurocs))) {
    if (all(is.na(aurocs[, i]))) { next }
    top_candidates <- find_top_candidate(votes[, i], aurocs[, i])
    result[top_candidates$best, i] <- top_candidates$score
    result[top_candidates$second, i] <- top_candidates$second_score
    result[top_candidates$third, i] <- top_candidates$third_score
    result[top_candidates$forth, i] <- top_candidates$forth_score
    result[top_candidates$fifth, i] <- top_candidates$fifth_score
      
#     result[top_candidates$sixth, i] <- top_candidates$sixth_score
#     result[top_candidates$seventh, i] <- top_candidates$seventh_score
#     result[top_candidates$eighth, i] <- top_candidates$eighth_score
#     result[top_candidates$ninth, i] <- top_candidates$ninth_score
#     result[top_candidates$tenth, i] <- top_candidates$tenth_score
    
  }
  return(result)
}

find_top_candidate <- function(votes, aurocs) {
  candidates <- extract_top_candidates(aurocs, 10)
  n_candidates <- length(candidates)
  
  # 初始化 AUROC 矩阵来存储每对候选者的 AUROC 值
  auroc_matrix <- matrix(NA, nrow = n_candidates, ncol = n_candidates,
                         dimnames = list(candidates, candidates))
  
  for (i in seq_len(n_candidates)) {
    for (j in seq(2, n_candidates)) {
      best <- candidates[i]
      contender <- candidates[j]
      votes_best <- votes[names(votes) == best]
      votes_contender <- votes[names(votes) == contender]
      
      # 检查名称是否匹配
#       if (length(votes_best) == 0 | length(votes_contender) == 0) {
#         next
#       }
      
      # 计算 AUROC 值
      auroc <- compute_aurocs(
        as.matrix(c(votes_best, votes_contender)),
        as.matrix(rep(c(1, 0), c(length(votes_best), length(votes_contender))))
      )
      
      # 更新 AUROC 矩阵
      auroc_matrix[best, contender] <- auroc
      auroc_matrix[contender, best] <- 1 - auroc  # 对称矩阵
    }
  }
  
  # 从 AUROC 矩阵中选择最佳候选者、次优候选者和第三候选者
  mean_aurocs <- rowMeans(auroc_matrix, na.rm = TRUE)
  
  # 对候选者进行排序
  sorted_indices <- order(-mean_aurocs)
  sorted_candidates <- names(mean_aurocs)[sorted_indices]
  
  # 提取最佳、次优和第三候选者
  best <- sorted_candidates[1]
  second_best <- sorted_candidates[2]
  third_best <- sorted_candidates[3]
  forth_best <- sorted_candidates[4]
  fifth_best <- sorted_candidates[5]
#   sixth_best <- sorted_candidates[6]
#   seventh_best <- sorted_candidates[7]
#   eighth_best <- sorted_candidates[8]
#   ninth_best <- sorted_candidates[9]
#   tenth_best <- sorted_candidates[10]  
    
    
  
  # 计算最终的 AUROC 值
  final_score <- mean_aurocs[best]
  second_score <- mean_aurocs[second_best]
  third_score <- mean_aurocs[third_best]
  forth_score <- mean_aurocs[forth_best]
  fifth_score <- mean_aurocs[fifth_best]

#   sixth_score <- mean_aurocs[sixth_best]
#   seventh_score <- mean_aurocs[seventh_best]
#   eighth_score <- mean_aurocs[eighth_best]
#   ninth_score <- mean_aurocs[ninth_best]
#   tenth_score <- mean_aurocs[tenth_best]
    
    
  
  return(list(score = final_score, best = best, second = second_best, second_score = second_score,
              third = third_best, third_score = third_score,
             forth =forth_best , forth_score =forth_score,
             fifth= fifth_best, fifth_score = fifth_score #,
#              sixth = sixth_best,sixth_score=sixth_score,
#              seventh = seventh_best,seventh_score = seventh_score,
#              eighth = eighth_best,eighth_score = eighth_score,
#              ninth = ninth_best,ninth_score = ninth_score,
#              tenth =  tenth_best,tenth_score = tenth_score
             ))
}

extract_top_candidates <- function(aurocs, n = 10) {
  return(names(utils::head(sort(aurocs, decreasing=TRUE), n = n)))
}


#' Runs unsupervised version of MetaNeighbor
#'
#' When it is difficult to know how cell type labels compare across datasets this
#' function helps users to make an educated guess about the overlaps without
#' requiring in-depth knowledge of marker genes
#'
#' @param var_genes vector of high variance genes.
#' @param dat SummarizedExperiment object containing gene-by-sample expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#' @param trained_model default value NULL; a matrix containing a trained model
#' generated from MetaNeighbor::trainModel. If not NULL, the trained model is
#' treated as training data and dat is treated as testing data. If a trained model
#' is provided, fast_version will automatically be set to TRUE and var_genes will
#' be overridden with genes used to generate the trained_model
#' @param fast_version default value FALSE; a boolean flag indicating whether
#' to use the fast and low memory version of MetaNeighbor
#' @param node_degree_normalization default value TRUE; a boolean flag indicating
#' whether to use normalize votes by dividing through total node degree.
#' @param one_vs_best default value FALSE; a boolean flag indicating whether
#' to compute AUROCs based on a best match against second best match setting
#' (default version is one-vs-rest). This option is currently only relevant when fast_version = TRUE.
#' @param symmetric_output default value TRUE; a boolean flag indicating whether to average AUROCs in the output matrix.
#'
#' @return The output is a cell type-by-cell type mean AUROC matrix, which is built by treating each pair of cell types as testing and training data for
#' MetaNeighbor, then taking the average AUROC for each pair (NB scores will not
#' be identical because each test cell type is scored out of its own dataset,
#' and the differential heterogeneity of datasets will influence scores).
#' If symmetric_output is set to FALSE, the training cell types are displayed as columns and the test cell types are displayed as rows.
#' If trained_model was provided, the output will be a cell type-by-cell
#' type AUROC matrix with training cell types as columns and test cell types
#' as rows (no swapping of test and train, no averaging).
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes,
#'                              dat = mn_data,
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' celltype_NV
#'

#' @export
MetaNeighborUS <- function(var_genes = c(), dat, i = 1, study_id, cell_type,
                           trained_model = NULL, fast_version = FALSE,
                           node_degree_normalization = TRUE, one_vs_best = FALSE,
                           symmetric_output = TRUE) {

    dat    <- SummarizedExperiment::assay(dat, i = i)
    samples <- colnames(dat)
    if (!is.null(trained_model)) {
        trained_model <- as.matrix(trained_model)
        var_genes <- rownames(trained_model)[-1]
    }
    
    #check obj contains study_id
    if(length(study_id)!=length(samples)){
        stop('study_id length does not match number of samples')
    }

    #check obj contains cell_type
    if(length(cell_type)!=length(samples)){
        stop('cell_type length does not match number of samples')
    }

    matching_vargenes <- match(rownames(dat), var_genes)
    matching_vargenes_count   <- sum(!is.na(matching_vargenes))

    if(matching_vargenes_count < 2){
        stop("matching_vargenes should have more than 1 matching genes!",
             call. = TRUE)
    } else if(matching_vargenes_count < 5) {
        warning("matching_vargenes should have more matching genes!",
                immediate. = TRUE)
    }
    dat <- dat[!is.na(matching_vargenes),]
    if (!is.null(trained_model)) {
        trained_model <- trained_model[c("n_cells", rownames(dat)),]
    }

    study_id <- as.character(study_id)
    cell_type <- as.character(cell_type)

    if (is.null(trained_model)) {
        if (fast_version) {
          cell_NV <- MetaNeighborUSLowMem(dat, study_id, cell_type,
                                          node_degree_normalization, one_vs_best)
        } else {
          cell_NV <- MetaNeighborUSDefault(dat, study_id, cell_type,
                                           node_degree_normalization)
        }
        if (symmetric_output) {
            cell_NV <- (cell_NV+t(cell_NV))/2
        }
    } else {
        cell_NV <-  MetaNeighborUS_from_trained(trained_model, dat, study_id, cell_type,
                                                node_degree_normalization, one_vs_best)
    }
    return(cell_NV)
}

MetaNeighborUSDefault <- function(dat, study_id, cell_type, node_degree_normalization = TRUE) {
    dat <- as.matrix(dat)
    pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
    pheno$StudyID_CT <- makeClusterName(pheno$study_id, pheno$cell_type)
    celltypes   <- unique(pheno$StudyID_CT)
    cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
    rownames(cell_labels) <-colnames(dat)
    colnames(cell_labels) <- celltypes

    for(i in seq_along(celltypes)){
        type <- celltypes[i]
        matching_celltype <- match(pheno$StudyID_CT, type)
        cell_labels[!is.na(matching_celltype),i]  <- 1
    }

    cor_data    <- stats::cor(dat, method="s")
    rank_data   <- cor_data*0
    rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
    rank_data[is.na(rank_data)] <- 0
    rank_data   <- rank_data/max(rank_data)
    sum_in      <- (rank_data) %*% cell_labels
    
    if (node_degree_normalization) {
        sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum),
                              ncol = dim(sum_in)[2],
                              nrow = dim(sum_in)[1])
        predicts    <- sum_in/sum_all
    } else {
        predicts <- sum_in        
    }

    cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
    colnames(cell_NV) <- colnames(cell_labels)
    rownames(cell_NV) <- colnames(cell_labels)

    for(i in seq_len(dim(cell_labels)[2])){
        predicts_temp <- predicts

        matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
        unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
        matching_studyID  <- match(pheno$study_id, unique_studyID)
        pheno2            <- pheno[!is.na(matching_studyID),]
        predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
        predicts_temp     <- apply(abs(predicts_temp),
                                   MARGIN = 2,
                                   FUN = rank,
                                   na.last= "keep",
                                   ties.method="average")


        filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
        matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
        filter[!is.na(matches),seq_along(celltypes)] <- 1

        negatives = which(filter == 0, arr.ind = TRUE)
        positives = which(filter == 1, arr.ind = TRUE)

        predicts_temp[negatives] <- 0

        np <- colSums(filter, na.rm = TRUE)
        nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
        p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)

        cell_NV[i,]= (p/np - (np+1)/2)/nn
    }
    return(cell_NV)
}

# The fast version is vectorized according to the following equations
# (Note that the point of these equations is to *never* compute the cell-cell network
#  by reordering the matrix operations):
#  - INPUTS:
#    + Q = test (Query) data (genes x cells)
#    + R = train (Ref) data (genes x cells)
#    + L = binary encoding of training cell types (Labels) (cells x cell types)
#    + S = binary encoding of train Studies (cells x studies)
#  - NOTATIONS:
#    + X* = normalize_cols(X) ~ scale(colRanks(X)) denotes normalized data
#           (Spearman correlation becomes a simple dot product on normalized data)
#    + N = Spearman(Q,R) = t(Q*).R* is the cell-cell similarity network
#    + CL = R*.L are the cell type centroids (in the normalized space)
#    + CS = R*.S are the study centroids (in the normalized space)
#    + 1.L = colSums(L) = number of cells per (train) cell type
#    + 1.S = colSums(S) = number of cells per (train) study
#  - WITHOUT node degree normalization
#    + Votes = N.L = t(Q*).R*.L = t(Q*).CL
#  - WITH node degree normalization
#    + Network becomes N+1 to avoid negative values
#    + Votes = (N+1).L = N.L + 1.L = t(Q*).CL + 1.L
#    + Node degree = (N+1).S = t(Q*).CS + 1.S
#    + Note: Node degree is computed independently for each train study.
#
MetaNeighborUSLowMem <- function(dat, study_id, cell_type,
                                 node_degree_normalization = TRUE, one_vs_best = FALSE) {
  dat <- normalize_cols(dat)
  label_matrix <- design_matrix(makeClusterName(study_id, cell_type))
  is_na <- matrixStats::colAnyNAs(dat)
  dat <- dat[, !is_na]
  label_matrix <- label_matrix[!is_na,]
  cluster_centroids <- dat %*% label_matrix
  n_cells_per_cluster <- colSums(label_matrix)
  result <- predict_and_score(dat, study_id[!is_na], cell_type[!is_na],
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization, one_vs_best)
  result <- result[, rownames(result)]
  return(result)
}

# 
                    
                    
# WARNING: function assumes that data have been normalized with normalize_cols
predict_and_score <- function(dat, study_id, cell_type,
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization = TRUE, one_vs_best = FALSE) {
  colnames(dat) <- makeClusterName(study_id, cell_type)
  votes <- crossprod(dat, cluster_centroids)
  if (node_degree_normalization) {
    votes <- normalize_node_degree(votes, dat, cluster_centroids, n_cells_per_cluster)
  }
  result <- c()
  for (test_study in unique(study_id)) {
    study_votes <- votes[study_id == test_study,]
    aurocs <- compute_aurocs(study_votes, design_matrix(rownames(study_votes)))
    if (one_vs_best) {
      result <- rbind(result, compute_1v1_aurocs(study_votes, aurocs)) # 我理解的这个输入的就是aurocs的值，但是是按照每行进行计算，因为是rbind垂直堆叠
    } else {
      result <- rbind(result, aurocs)
    }
  }
  return(result)
}

normalize_node_degree <- function(votes, dat, cluster_centroids, n_cells_per_cluster) {
    # node degree is normalized by train study -> compute study centroids
    centroid_study_label <- getStudyId(colnames(cluster_centroids))
    study_matrix <- design_matrix(centroid_study_label)
    study_centroids <- cluster_centroids %*% study_matrix
    n_cells_per_study <- n_cells_per_cluster %*% study_matrix
    train_study_id <- colnames(study_matrix)
    
    node_degree <- crossprod(dat, study_centroids)
    
    # shift to positive values and normalize node degree
    result <- sweep(votes, 2, n_cells_per_cluster, "+")
    node_degree <- sweep(node_degree, 2, n_cells_per_study, "+")
    for (train_study in unique(train_study_id)) {
      is_train <- centroid_study_label == train_study
      result[, is_train] <- result[, is_train] / node_degree[, train_study]
    }
  return(result)
}
                    
MetaNeighborUS_from_trained <- function(trained_model, test_dat, study_id, cell_type,
                                        node_degree_normalization = TRUE, one_vs_best = FALSE) {
  dat <- normalize_cols(test_damakeClusterNamet)
  is_na <- matrixStats::colAnyNAs(dat)
  dat <- dat[, !is_na]
  cluster_centroids <- trained_model[-1,]
  n_cells_per_cluster <- trained_model[1,]
    
  result <- predict_and_score(dat, study_id[!is_na], cell_type[!is_na],
                              cluster_centroids, n_cells_per_cluster,
                              node_degree_normalization, one_vs_best)
  return(result)
}
