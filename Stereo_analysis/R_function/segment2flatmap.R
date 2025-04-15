Archr_col <- c('#D51F26','#272E6A','#208A42','#89288F','#F47D2B','#FEE500','#8A9FD1',
               '#C06CAB','#E6C2DC','#90D5E4','#89C75F','#F37B7D','#9983BD','#D24B27',
               '#3BBCA8','#6E4B9E','#0C727C','#7E1416','#D8A767')
caculate_length <- function(tmp){
    colnames(tmp) <- c("x", "y")
    tmp1 <- tmp[1:(nrow(tmp) - 1), c("x", "y")]
    tmp1[, c("x1", "y1")] <- tmp[2:nrow(tmp), c("x", "y")]
    return(sum(((tmp1$x - tmp1$x1)^2 + (tmp1$y - tmp1$y1)^2)^0.5))}

trakem_transform <- function(x,y,matrix,matrix_bin,transform_bin){
    scale = matrix(c(transform_bin/matrix_bin,0,0,0,transform_bin/matrix_bin,0,0,0,1),ncol=3,nrow=3)
    matrix_t <- rbind(matrix(c(as.numeric(strsplit(matrix,',')[[1]])),nrow=2),c(0,0,1))
    matrix_t <- solve(scale) %*% matrix_t %*% scale
    coor_matrix <- cbind(x,y,1)
    coor_matrix_t <- t(matrix_t %*% t(coor_matrix))
    rx <- coor_matrix_t[,1]
    ry <- coor_matrix_t[,2]
    return(list('rx'=rx, 'ry'=ry))
    }
creat_segment_obj <- function(input_data,cluster1=c('Gaba','Glut')){
    input_data <- input_data[input_data$cluster1%in%cluster1,]
    input_data$label <- paste0(input_data$slice,'_segment_',input_data$segment)
    input_data$geneID <- input_data$celltype
    input_data$MIDCount <- 1
    
    gene=1:length(unique(input_data$geneID))
    names(gene)=unique(input_data$geneID)
    cell=1:length(unique(input_data$label))
    names(cell)=unique(input_data$label)
    mat1=sparseMatrix(i = gene[input_data$geneID],j=cell[ as.character(input_data$label) ], x=input_data$MIDCount)
    rownames(mat1)=names(gene)
    colnames(mat1)=names(cell)
    obj2=CreateSeuratObject(counts = mat1)
    return(obj2)
}
prepare_obj <- function(tmp){
    tmp <- NormalizeData(tmp,verbose = F)
    tmp@assays$RNA@var.features <- rownames(tmp)  
    tmp <- ScaleData(tmp,verbose = F)
    tmp <- RunPCA(tmp,verbose = F,npcs=10)
    return(tmp)
}
imputation_celltype <- function(obj,segment.meta.data,assays='Neuron'){
    #imputation
    obj@meta.data <- obj@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA')]
    #create zeros matrix which segment is not in segment meta data
    imputation_segment <- segment.meta.data$id[!segment.meta.data$id %in% colnames(obj)]
    imputation_matrix <- matrix(data=0,ncol=length(imputation_segment),nrow=length(rownames(obj)))
    colnames(imputation_matrix) <- imputation_segment
    rownames(imputation_matrix) <- rownames(obj)
    obj_imputation <- CreateSeuratObject(imputation_matrix)
    obj_imputation$imputation <- 'imputation'
    obj$imputation <- 'Raw'
    #combined raw matrix and imputation matrix
    obj_all <- merge(obj,obj_imputation)
    
    rownames(segment.meta.data) <- segment.meta.data$id
    obj_all@meta.data[,c('id','slice','segment','X_axis','Y_axis','Y_axis_smooth','length','region','VL_region')] <-  
    segment.meta.data[colnames(obj_all),c('id','slice','segment','X_axis','Y_axis','Y_axis_smooth','length','region_imputation','VL_smooth')]
    obj_all <- NormalizeData(obj_all)
    obj <- obj_all
    #imputation each celltypes by nearest 9 point from raw matrix
    for( i in rownames(obj)){
        obj@assays$RNA@data[i,obj$imputation=='imputation'] <- smooth_kNN(
            obj@meta.data[obj$imputation=='Raw',c('X_axis','Y_axis_smooth')],
            obj@meta.data[obj$imputation=='imputation',c('X_axis','Y_axis_smooth')],
            obj@assays$RNA@data[i,obj$imputation=='Raw'],
            knn=9,round=1)
    }
    obj@assays$RNA@var.features <- rownames(obj)
    obj[[assays]] <- obj[['RNA']]
    obj@active.assay <- assays
    obj[['RNA']] <- NULL
    
    return(obj)
}

smooth_celltype <- function(obj,assays='Neuron'){
    #smooth data use smooth_kNN
    obj[[paste0(assays,'_','smooth')]] <- obj[[assays]]
    data_matrix <- as.matrix(obj[[paste0(assays,'_','smooth')]]@data)
    for( i in rownames(obj)){
        data_matrix[i,] <- smooth_kNN(
            obj@meta.data[,c('X_axis','Y_axis_smooth')],
            obj@meta.data[,c('X_axis','Y_axis_smooth')],
            data_matrix[i,],
            knn=25,round=10)
            }
        data_matrix <- CreateAssayObject(data=data_matrix)
        data_matrix <- CreateSeuratObject(data_matrix)
        obj[[paste0(assays,'_','smooth')]]@data <- data_matrix@assays$RNA@data
        obj@active.assay <- paste0(assays,'_','smooth')
        obj <- ScaleData(obj)
        return(obj)
}