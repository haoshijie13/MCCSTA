#Visualize function
library(ggplot2)
library(grid)
Archr_col <- c('#D51F26','#272E6A','#208A42','#89288F','#F47D2B','#FEE500','#8A9FD1',
               '#C06CAB','#E6C2DC','#90D5E4','#89C75F','#F37B7D','#9983BD','#D24B27',
               '#3BBCA8','#6E4B9E','#0C727C','#7E1416','#D8A767')
plot_col_exp <- function(obj,col_exp,rev=FALSE,vmid=FALSE,min.cutoff=0,border=TRUE,max.cutoff=1,color=RColorBrewer::brewer.pal(11,'Spectral')[-6]){
    exp <- obj[,col_exp]
    if(min.cutoff!=0 | max.cutoff!=1){
    max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
    min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
    exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
    exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    color_list <- rev(color)
    if(rev){color_list <- color}
    p1 <- ggplot()+
    geom_tile(data=obj,aes(y=Y_axis,x=X_axis,fill=exp,height=length))+
    scale_fill_gradientn(colours = color_list,name=col_exp,na.value = 'gray')+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
    theme_void()+
    coord_fixed()
    if(vmid!=FALSE){
        p1 <- ggplot()+
        geom_tile(data=obj,aes(y=Y_axis,x=X_axis,fill=exp,height=length+0.1))+
        scale_fill_gradientn(colours = color_list,name=col_exp,
                            values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))))+
        theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
        theme_void()+
        coord_fixed()}
    if(border){p1 <- p1+geom_segment(data=Border,aes(x=X,y=Y,xend=X1,yend=Y1),lwd=0.8)}
    return(p1)
    }
plot_features <- function(obj,features,min.cutoff=0,border=TRUE,max.cutoff=1,assays='smooth',slot='scale.data',vmid=0.5,color=c('navy','red','yellow')){
plot.list <- lapply(features,function(i){
    exp <- slot(obj@assays[[assays]],slot)[i,]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    p1 <- ggplot()+
    geom_tile(data=obj@meta.data,aes(y=Y_axis,x=X_axis,fill=exp,height=length))+
    scale_fill_gradientn(colours = color,name=i,na.value = 'gray',
                         values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))))+
    theme_void()+
    coord_fixed()
    #geom_segment(data=Border,aes(x=X,y=Y,xend=X1,yend=Y1),color='black',lwd=0.8)
    if(border){p1 <- p1+geom_segment(data=Border,aes(x=X,y=Y,xend=X1,yend=Y1),color='black',lwd=0.8)}
    return(p1)
        })
    p <- plot.list[[1]]
    if(length(features)>1){
    for(i in 2:length(features)){p <- p+plot.list[[i]]}
    return(p)
        }else{
    return(p)
        }

}
save_png_plot <- function(plot,dir_path,name,height=200,width=500,save='png'){
    if(!dir.exists(dir_path)){dir.create(dir_path)}
    if(save=='ggsave'){ggsave(paste0(dir_path, "/", name, ".png"),plot + NoLegend(), height = height, width = width,bg = "transparent", units = 'px')
                      }else if(save=='png'){
        png(paste0(dir_path, "/", name, ".png"), height = height, width = width,bg = "transparent")
        print(plot + NoLegend())
        dev.off()}
    legend <- cowplot::get_legend(plot)
    pdf(paste0(dir_path,'/',name,'.legend.pdf'))
    print(grid::grid.draw(legend))
    dev.off()
}