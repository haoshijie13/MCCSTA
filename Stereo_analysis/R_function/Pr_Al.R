#calculate segment correlation with choice region, cal_matrix(segment × features), region_matrix(region × features)
region_cor_segment_one <- function(cal_matrix,region_matrix,choice_region){
    use_row <- intersect(rownames(cal_matrix),rownames(region_matrix))
    tmp_cor <- cor(cal_matrix[use_row,, drop = FALSE],region_matrix[use_row,choice_region, drop = FALSE])
    if(ncol(tmp_cor)>1){tmp_cor <- matrixStats::rowMaxs(tmp_cor)}else{tmp_cor <- tmp_cor[,choice_region]}
    tmp_cor <- Mean_normalize(tmp_cor)/2
    return(list('cor_value'=tmp_cor))
}
#calculate Pr-Al Index by input Pr_matrix (Pr × features) or Al_matrix (Al × features)
PrAl_index_cal <- function(cal_matrix,Pr_matrix,Al_matrix,normalize=FALSE){
    use_row <- intersect(rownames(Pr_matrix),rownames(Al_matrix))
    use_row <- intersect(use_row,rownames(cal_matrix))
    tmp_cor <- cor(cal_matrix[use_row,, drop = FALSE],cbind(Pr_matrix[use_row,, drop = FALSE],Al_matrix[use_row,, drop = FALSE]))
    if(ncol(Pr_matrix)>1){Cor_Pr <- matrixStats::rowMaxs(tmp_cor[,colnames(Pr_matrix)])}else{Cor_Pr <- tmp_cor[,colnames(Pr_matrix)]}
    if(ncol(Al_matrix)>1){Cor_Al <- matrixStats::rowMaxs(tmp_cor[,colnames(Al_matrix)])}else{Cor_Al <- tmp_cor[,colnames(Al_matrix)]}
    Cor_Pr <- (Cor_Pr - min(Cor_Pr))/(max(Cor_Pr) - min(Cor_Pr))
    Cor_Al <- (Cor_Al - min(Cor_Al))/(max(Cor_Al) - min(Cor_Al))
    if(normalize){
    Cor_Pr <- Mean_normalize(Cor_Pr)/2
    Cor_Al <- Mean_normalize(Cor_Al)/2}
    Pr_Al_index <- (Cor_Pr-Cor_Al)/((1-Cor_Pr)+(1-Cor_Al))
    return(list('Pr-Al_index'=Pr_Al_index,'Cor_Pr'=Cor_Pr,'Cor_Al'=Cor_Al))
}

prepare_obj <- function(tmp){
    tmp <- NormalizeData(tmp,verbose = F)
    if('layer_enrich'%in%colnames(tmp@assays$Neuron@meta.features)){
        tmp@assays$Neuron@var.features <- rownames(tmp)[tmp@assays$Neuron@meta.features$layer_enrich]
    }else(tmp@assays$Neuron@var.features <- rownames(tmp))
    tmp <- ScaleData(tmp,verbose = F)
    tmp <- RunPCA(tmp,verbose = F,npcs=5)
    tmp@meta.data[,c('Cell-PC1')] <- tmp@reductions$pca@cell.embeddings[,1]
    return(tmp)}

prepare_obj_gene <- function(obj,species){
    tmp <- obj
    tmp <- CreateSeuratObject(obj[['RNA']]@counts[rownames(obj)[obj[['RNA']]@meta.features$homo_exp],])
    tmp <- NormalizeData(tmp,verbose = F)
    tmp@assays$RNA@var.features <- rownames(tmp)
    tmp <- ScaleData(tmp,verbose = F)
    tmp <- RunPCA(tmp,verbose = F,npcs=5)
    obj@meta.data[colnames(tmp),c('Gene-PC1')] <- tmp@reductions$pca@cell.embeddings[,1]
    obj[['pca']] <- tmp[['pca']]
    obj <- NormalizeData(obj,verbose = F)
    return(obj)}

#calculate correlation between celltype distribution with Pr-Al Index  
CI_cal <- function(obj){
    cell_cor <- sapply(rownames(obj@assays$Neuron),function(i){return(cor(obj$`Pr-Al-Index`,obj@assays$Neuron@scale.data[i,]))})
    return(cell_cor)}
CI_cal_p <- function(obj){
    cell_p <- sapply(rownames(obj@assays$Neuron),function(i){return(cor.test(obj$`Pr-Al-Index`,obj@assays$Neuron@scale.data[i,])$p.value)})
    return(cell_p)}

#covert Pr-Al Index to Intersection Index
intersection_convert <- function(input_vector_raw){
    input_vector <- input_vector_raw
    use_thre <- 0
    up_vector <- input_vector[input_vector<use_thre]
    down_vector <- input_vector[input_vector>=use_thre]
    input_vector[input_vector_raw<use_thre] <- ((up_vector-min(up_vector))/(max(up_vector)-min(up_vector)))
    input_vector[input_vector_raw>=use_thre] <- 1-(down_vector-min(down_vector))/(max(down_vector)-min(down_vector))
    return(input_vector)
}

plot_cor_raster <- function(obj_segment,X,Y,var,fitting=FALSE){
    p <- ggplot(data=obj_segment@meta.data,aes(x=obj_segment@meta.data[,X],y=obj_segment@meta.data[,Y]))+
    scattermore::geom_scattermore(aes(color=obj_segment@meta.data[,var]),
                     pointsize = 2)+
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")[-6]),name=var)+
    xlab(X)+
    ylab(Y)+
    theme_classic()+
    theme(aspect.ratio=1,text=element_text(size=16))+
    ggtitle(round(cor(obj_segment@meta.data[,X],obj_segment@meta.data[,Y]),2))
    if(fitting){p <- p+geom_smooth(method = 'lm',lwd=2,se = T,color='gray')}
    return(p)
}

plot_cor_dot <- function(obj,X,Y,var,fitting=FALSE,method='pearson'){
    p <- ggplot(data=obj@meta.data,aes(x=obj@meta.data[,X],y=obj@meta.data[,Y]))+
    geom_point(aes(color=obj@meta.data[,var]),show.legend = T, shape=21, size = 3, stroke = 1.5)+
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")[-6]),name=var)+
    xlab(X)+
    ylab(Y)+
    theme_classic()+
    theme(aspect.ratio=1,text=element_text(size=16))+
    ggtitle(paste0('R=',round(cor(obj@meta.data[,X],obj@meta.data[,Y],method=method),2),' ',
                   'p=',signif(cor.test(obj@meta.data[,X],obj@meta.data[,Y],method=method)$p.value,3)))
    if(fitting){p <- p+geom_smooth(method = 'lm',lwd=2,se = T,color='gray')}
    return(p)
}