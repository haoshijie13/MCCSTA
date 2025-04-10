library(RANN)
#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#round is the average round
#knn is the n nearest point from reference
smooth_kNN <- function(initial_df,query_df,sm_vector,round=100,knn=30){
    result <- nn2(initial_df,query_df,k=knn)
    for(i in 1:round){
        if(i==1){initail_col <- sm_vector}else{initail_col <- sm_vector_sm}
        sm_vector_sm <- rowMeans(matrix(initail_col[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE))
    }
    return(sm_vector_sm)
}
#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#knn is the n nearest point from reference
#If the value counts exceed the threshold, specific actions apply; otherwise, use the nearest point as the query value.
winner_kNN <- function(initial_df,query_df,sm_vector,knn=5,threhold=1){
    result <- nn2(initial_df,query_df,k=knn)
    initail_col <- sm_vector
    sm_vector_sm <- matrix(initail_col[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE)
    sm_vector_sm <- sapply(1:nrow(sm_vector_sm),function(x){
        tmp_table <- table(sm_vector_sm[x,])
        name <- names(sort(tmp_table,decreasing = T))[1]
        if(tmp_table[name]>threhold){return(name)}else{return(sm_vector_sm[x,1])}
    })
    return(sm_vector_sm)
}
ROC_kNN <- function (initial_df, vector, knn = 9) 
{
    result <- nn2(initial_df, initial_df, k = knn)
    vector_n <- matrix(vector[result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    x_n <- matrix(initial_df[, 1][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    y_n <- matrix(initial_df[, 2][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    vector_max <- matrixStats::rowMaxs(vector_n)
    vector_min <- matrixStats::rowMins(vector_n)
    gradient_top <- abs(vector_max - vector_min)
    gradient_bottom <- sapply(1:nrow(vector_n), function(i) {
        x_max <- x_n[i, ][vector_n[i, ] == vector_max[i]][1]
        y_max <- y_n[i, ][vector_n[i, ] == vector_max[i]][1]
        x_min <- x_n[i, ][vector_n[i, ] == vector_min[i]][1]
        y_min <- y_n[i, ][vector_n[i, ] == vector_min[i]][1]
        distance <- ((x_max - x_min)^2 + (y_max - y_min)^2)^0.5
        #distance <- as.matrix(dist(cbind(c(x_max,x_min),c(y_max,y_min))))
        #distance <- mean(distance[1:length(x_max),(length(x_max)+1):(length(x_max)+length(x_min))])
        return(distance)
    })
    result <- gradient_top/gradient_bottom
    result[is.na(result)] <- 0
    result[is.infinite(result)] <- 0
    result <- (result-min(result))/(max(result)-min(result))
    result[is.na(result)] <- 0
    return(result)
}
