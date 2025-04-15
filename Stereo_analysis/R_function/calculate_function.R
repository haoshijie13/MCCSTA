#Calculate function
matrix_groupby <- function(matrix_input,group.by,type='row',cal='sum'){
    if(type=='col'){matrix_input <- t(matrix_input)}
    combined_matrix <- data.frame('group'=group.by,'raw'=rownames(matrix_input),'identify'=1)
    combined_matrix <- reshape2::acast(combined_matrix,group~raw,fill = 0,value.var = 'identify')
    combined_matrix <- combined_matrix[,rownames(matrix_input)]
    if(cal=='sum'){matrix_output <- combined_matrix%*%matrix_input}else if(cal=='mean'){
        combined_matrix <- combined_matrix/rowSums(combined_matrix)
        matrix_output <- combined_matrix%*%matrix_input
    }
    if(type=='col'){matrix_output <- t(matrix_output)}
    return(matrix_output)}
#Calculate pearson C.C. after filter NA
cor_myself <- function(matrix_input,return.p=T){
    matrix_output <- matrix(0,nrow = ncol(matrix_input),ncol=ncol(matrix_input))
    rownames(matrix_output) <- colnames(matrix_input)
    colnames(matrix_output) <- colnames(matrix_input)
    if(return.p){matrix.p <- matrix_output}
    for(i in colnames(matrix_input)){
        for(j in colnames(matrix_input)){
            use_index <- !(is.na(matrix_input[,i])|is.na(matrix_input[,j]))
            result <- cor.test(matrix_input[use_index,i],matrix_input[use_index,j])
            matrix_output[i,j] <- result$estimate
            if(return.p){matrix.p[i,j] <- result$p.value}
            }
    }
    if(return.p){return(list('R'=matrix_output,'p.value'=matrix.p))}else(return(matrix_output))
}
cor_twomatrix <- function(matrix_1,matrix_2,return.p=F){
    matrix_output <- matrix(0,nrow = ncol(matrix_1),ncol=ncol(matrix_2))
    rownames(matrix_output) <- colnames(matrix_1)
    colnames(matrix_output) <- colnames(matrix_2)
    matrix.p <- matrix_output
    n <- 1
    all_n <- ncol(matrix_1)*ncol(matrix_2)
    pb <- txtProgressBar(style=3)
    for(i in colnames(matrix_1)){
        for(j in colnames(matrix_2)){
            use_index <- !(is.na(matrix_1[,i])|is.na(matrix_2[,j]))
            result <- cor.test(matrix_1[use_index,i],matrix_2[use_index,j])
            matrix_output[i,j] <- result$estimate
            if(return.p){matrix.p[i,j] <- result$p.value}
            setTxtProgressBar(pb, n/all_n)
            n <- n+1
            }
    }
    close(pb)
    if(return.p){return(list('R'=matrix_output,'p.value'=matrix.p))}else(return(matrix_output))
}
Mean_normalize <- function(input_vector_raw){
    input_vector <- input_vector_raw
    up_vector <- input_vector[input_vector_raw>=mean(input_vector_raw)]
    down_vector <- input_vector[input_vector_raw<mean(input_vector_raw)]
    input_vector[input_vector_raw>=mean(input_vector_raw)] <- ((up_vector-min(up_vector))/(max(up_vector)-min(up_vector)))+1
    input_vector[input_vector_raw<mean(input_vector_raw)] <- (down_vector-min(down_vector))/(max(down_vector)-min(down_vector))
    return(input_vector)
}
normalize_columns <- function(mat) {
  min_vals <- apply(mat, 2, min)
  max_vals <- apply(mat, 2, max)
  normalized_mat <- t((t(mat) - min_vals) / (max_vals - min_vals))
  return(normalized_mat)
}
normalize_row <- function(mat) {
  min_vals <- apply(mat, 1, min)
  max_vals <- apply(mat, 1, max)
  normalized_mat <- (mat - min_vals) / (max_vals - min_vals)
  return(normalized_mat)
}
order_matrix <- function(input_matrix,row_order){
    input_matrix <- scale(input_matrix)
    input_matrix[is.na(input_matrix)] <- 0
    return_df <- dplyr::bind_rows(lapply(colnames(input_matrix),function(x){
        max_value <- max(input_matrix[,x])
        max_label <- rownames(input_matrix)[input_matrix[,x]==max_value][1]
        df <- data.frame(var=x,label=max_label,value=max_value)
        return(df)}))
    rownames(return_df) <- return_df$var
    return_df <- dplyr::bind_rows(lapply(row_order,function(x){
        tmp_df <- return_df[return_df$label==x,]
        tmp_df <- tmp_df[order(tmp_df$value,decreasing =T),]
        return(tmp_df)
    }))
    return(return_df)
}