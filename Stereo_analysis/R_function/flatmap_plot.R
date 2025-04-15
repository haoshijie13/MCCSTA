library(ggplot2)
library(sf)
Sys.setenv(PROJ_LIB = "/mnt/SSD4Ta/home/huangzhi/anaconda3/envs/work/share/proj")

flat_mar <- read_sf('data/MRI_shp_file/Marmoset.shp')
flat_mac <- read_sf('data/MRI_shp_file/Macaque.shp')
flat_hum <- read_sf('data/MRI_shp_file/Human.shp')
flat_mou <- read_sf('data/MRI_shp_file/Mouse.shp')

plot_mar_flatmap_col <- function(obj.meta,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    p1 <- ggplot(flat_mar) +
    geom_sf(aes(fill=obj.meta[flat_mar$region,plot_names]),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}
plot_mac_flatmap_col <- function(obj.meta,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    tmp <- obj.meta
    tmp[tmp$region%in%c('Iam','Iapm'),'region'] <- 'Iam-Iapm'
    tmp[tmp$region%in%c('PEc','PEci'),'region'] <- 'PEc-PEci'
    tmp[tmp$region%in%c('Ia','Id'),'region'] <- 'Ia-Id'
    tmp[tmp$region%in%c('PG','Opt','DP'),'region'] <- 'PG-Opt-DP'
    tmp[tmp$region%in%c('36r','36p'),'region'] <- '36r-36p'
    tmp[tmp$region%in%c('V4','V4t'),'region'] <- 'V4-V4t'
    tmp[tmp$region%in%c('PF','PFG'),'region'] <- 'PF-PFG'
    tmp[tmp$region%in%c('13a','13b'),'region'] <- '13a-13b'
    tmp[tmp$region%in%c('EO','EI','ELr','ELc'),'region'] <- 'EO-EI-ELr-ELc'
    tmp <- aggregate(list('feature'=tmp[,plot_names]),by=list('region'=tmp$region),mean)
    rownames(tmp) <- tmp$region
    
    p1 <- ggplot(flat_mac) +
    geom_sf(aes(fill=tmp[flat_mac$region,]$feature),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_hum_flatmap_col <- function(obj.meta,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    p1 <- ggplot(flat_hum) +
    geom_sf(aes(fill=obj.meta[flat_hum$region,plot_names]),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_mou_flatmap_col <- function(obj.meta,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    p1 <- ggplot(flat_mou) +
    geom_sf(aes(fill=obj.meta[flat_mou$region,plot_names]),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
return(p1)
}

#plot_features
plot_mar_flatmap_feature <- function(obj,feature,assays='RNA',slot='scale.data',vmid=0.5,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){ 
    obj$exp <- slot(obj@assays[[assays]],slot)[feature,]
    p1 <- ggplot(flatmap_rotation(flat_mar,-45)) +
    geom_sf(aes(fill=obj@meta.data[flat_mar$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}
plot_mac_flatmap_feature <- function(obj,feature,assays='RNA',slot='scale.data',vmid=0.5,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    obj$exp <- slot(obj@assays[[assays]],slot)[feature,]
    tmp <- obj@meta.data
    tmp[tmp$region%in%c('Iam','Iapm'),'region'] <- 'Iam-Iapm'
    tmp[tmp$region%in%c('PEc','PEci'),'region'] <- 'PEc-PEci'
    tmp[tmp$region%in%c('Ia','Id'),'region'] <- 'Ia-Id'
    tmp[tmp$region%in%c('PG','Opt','DP'),'region'] <- 'PG-Opt-DP'
    tmp[tmp$region%in%c('36r','36p'),'region'] <- '36r-36p'
    tmp[tmp$region%in%c('V4','V4t'),'region'] <- 'V4-V4t'
    tmp[tmp$region%in%c('PF','PFG'),'region'] <- 'PF-PFG'
    tmp[tmp$region%in%c('13a','13b'),'region'] <- '13a-13b'
    tmp[tmp$region%in%c('EO','EI','ELr','ELc'),'region'] <- 'EO-EI-ELr-ELc'
    tmp <- aggregate(list('feature'=tmp[,'exp']),by=list('region'=tmp$region),mean)
    rownames(tmp) <- tmp$region
    
    p1 <- ggplot(flat_mac) +
    geom_sf(aes(fill=tmp[flat_mac$region,]$Gene_EI),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_hum_flatmap_feature <- function(obj,feature,assays='RNA',slot='scale.data',vmid=0.5,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    obj$exp <- slot(obj@assays[[assays]],slot)[feature,]
    tmp <- obj@meta.data
    p1 <- ggplot(flat_hum) +
    geom_sf(aes(fill=tmp[flat_hum$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_mou_flatmap_feature <- function(obj,feature,assays='RNA',slot='scale.data',vmid=0.5,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    obj$exp <- slot(obj@assays[[assays]],slot)[feature,]
    tmp <- obj@meta.data
    p1 <- ggplot(flat_mou) +
    geom_sf(aes(fill=tmp[flat_mou$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}
