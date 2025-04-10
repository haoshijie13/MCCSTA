cell_color <- c(colorRampPalette(c('#D51F26','red'))(2)[c(1,2)],
                colorRampPalette(c('#272E6A','navy'))(3)[c(1,2,3)],
                colorRampPalette(c('#208A42','#89C75F'))(5)[c(1,2,3,4,5)],
                colorRampPalette(c('#89288F','#C06CAB'))(4)[c(1,2,3,4)],
                colorRampPalette(c('#F47D2B','#D24B27'))(3)[c(1,2,3)],
                colorRampPalette(c('#FEE500','yellow'))(2)[c(1,2)],
                colorRampPalette(c('black','white'))(10)[c(4,6,8)]
                )
names(cell_color) <- c('VLMC','Ast',
                       'L2','RELN','VIP',
                       'L3','L2/3','L3/4/5','SST','LAMP5',
                       'L4','L3/4', 'PV-CHC','PV',        
                       'L5','L4/5','L5/6',
                       'L6','OLG',
                       'EC','MG','OPC')
Archr_col <- c('#D51F26','#272E6A','#208A42','#89288F','#F47D2B','#FEE500','#8A9FD1',
               '#C06CAB','#E6C2DC','#90D5E4','#89C75F','#F37B7D','#9983BD','#D24B27',
               '#3BBCA8','#6E4B9E','#0C727C','#7E1416','#D8A767')

plot_whole <- function(cell_data,slice,plot_function,dot_size=1.2,show_legend=FALSE,sag=FALSE){
    xmean <- (max(cell_data$rx) + min(cell_data$rx))/2
    ymean <- (max(cell_data$ry) + min(cell_data$ry))/2
    cell_data <- cell_data[cell_data$slice==slice,]
    cell_data <- cell_data[cell_data$cluster2!='discard',]
    
    
    if(plot_function=='cluster2'){
    cell_data <- rbind(cell_data[!cell_data$cluster2%in%c('L2','L3','L4','L5','L6'),],cell_data[cell_data$cluster2%in%c('L2','L3','L4','L5','L6'),])
    p1 <- ggplot()+
    geom_point(data=cell_data,
               aes(x=rx,y=ry,color=cluster2),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = cell_color[c('L2','L2/3','L3','L3/4','L3/4/5','L4','L4/5','L5','L5/6','L6',
                            'RELN','VIP','LAMP5','PV','PV-CHC','SST',
                                'VLMC','Ast','OLG','EC','MG','OPC')],
                      breaks=c('L2','L2/3','L3','L3/4','L3/4/5','L4','L4/5','L5','L5/6','L6',
                            'RELN','VIP','LAMP5','PV','PV-CHC','SST',
                                'VLMC','Ast','OLG','EC','MG','OPC'))+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))
    if(!sag){p1 <- p1+xlim(xmean-14000,xmean+14000)+ylim(ymean-19000,ymean+19000)
           }else{p1 <- p1+xlim(xmean-35000,xmean+35000)+ylim(ymean-20000,ymean+20000)}
    return(p1)
        
    }else if(plot_function=='layer')
    p1 <- ggplot()+
    geom_point(data=cell_data,
               aes(x=rx,y=ry,color=layer),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = c(Archr_col[1:6],Archr_col[8],Archr_col[7]),breaks=c('L1','L2','L3','L4','L5','L6','L2/3/4','WM'))+
    coord_fixed()+
    theme_void()+
    xlim(xmean-14000,xmean+14000)+
    ylim(ymean-19000,ymean+19000)+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))
    return(p1)
}
rotation_coor <- function(raw_df,rotation){
    r_pi <- rotation*pi/180
    r_matrix <- matrix(c(cos(r_pi),-sin(r_pi),
                         sin(r_pi),cos(r_pi)),ncol=2,nrow=2)
    r_df <- as.matrix(raw_df)%*%r_matrix
    return(r_df)
    }
plot_panel <- function(cell_data,slice,xmean,ymean,rotation,plot_function,celltype='None',color='red',dot_size=0.8,show_legend=FALSE){
    cell_data <- cell_data[cell_data$slice==slice,]
    cell_data[,c('rx','ry')] <- rotation_coor(cell_data[,c('rx','ry')],rotation)
    if(plot_function=='test'){
    p1 <- ggplot()+
    geom_point(data=cell_data[ seq(1,nrow(cell_data),length.out = 2000),],
               aes(x=rx,y=ry,color=layer),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = c(Archr_col[1:6],Archr_col[8],Archr_col[7]),breaks=c('L1','L2','L3','L4','L5','L6','L2/3/4','White'))+
    coord_fixed()+
    theme_classic()
    return(p1)       
                      }else if(plot_function=='celltype'){
    p1 <- ggplot()+
    geom_point(data=cell_data[!cell_data$celltype%in%celltype,],
               aes(x=rx,y=ry),color='gray',size=0.8,show.legend = show_legend)+
    geom_point(data=cell_data[cell_data$celltype%in%celltype,],
               aes(x=rx,y=ry),color=color,size=0.8,show.legend = show_legend)+
    #scale_color_lancet()+
    #scale_color_manual(values = pal_lancet())+
    xlim(xmean-4500,xmean+4500)+
    ylim(ymean-3000,ymean+3000)+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
    theme(legend.position = "bottom")
    return(p1)
        }else if(plot_function=='cluster2'){
    cell_data <- rbind(cell_data[!cell_data$cluster2%in%c('L2','L3','L4','L5','L6'),],cell_data[cell_data$cluster2%in%c('L2','L3','L4','L5','L6'),])
    p1 <- ggplot()+
    geom_point(data=cell_data,
               aes(x=rx,y=ry,color=cluster2),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = cell_color[c('L2','L2/3','L3','L3/4','L3/4/5','L4','L4/5','L5','L5/6','L6',
                            'RELN','VIP','LAMP5','PV','PV-CHC','SST',
                                'VLMC','Ast','OLG','EC','MG','OPC')],
                      breaks=c('L2','L2/3','L3','L3/4','L3/4/5','L4','L4/5','L5','L5/6','L6',
                            'RELN','VIP','LAMP5','PV','PV-CHC','SST',
                                'VLMC','Ast','OLG','EC','MG','OPC'))+
    xlim(xmean-4500,xmean+4500)+
    ylim(ymean-3000,ymean+3000)+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))
    return(p1)
        }else if(plot_function=='unsupervised_layer'){
    p1 <- ggplot()+
    geom_point(data=cell_data,
               aes(x=rx,y=ry,color=unsupervised_layer),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = c(Archr_col[1:6],Archr_col[8],Archr_col[7]),breaks=c('L1','L2','L3','L4','L5','L6','L2/3/4','White'))+
    xlim(xmean-4500,xmean+4500)+
    ylim(ymean-3000,ymean+3000)+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))
    return(p1)
        }else if(plot_function=='layer'){
    p1 <- ggplot()+
    geom_point(data=cell_data,
               aes(x=rx,y=ry,color=layer),size=dot_size,show.legend = show_legend)+
    scale_color_manual(values = c(Archr_col[1:6],Archr_col[8],Archr_col[7]),breaks=c('L1','L2','L3','L4','L5','L6','L2/3/4','White'))+
    xlim(xmean-4500,xmean+4500)+
    ylim(ymean-3000,ymean+3000)+
    coord_fixed()+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))
    return(p1)
        }
}
layer_color <-  c(Archr_col[1:6],Archr_col[8],Archr_col[7])
names(layer_color) <- c('L1','L2','L3','L4','L5','L6','L2/3/4','WM')
plot_layer_compare <- function(manual_label){
    manual_label_tmp <- manual_label[!manual_label$layer%in%c('WM','L1'),]
    Sample_plot_all <- dplyr::bind_rows(lapply(unique(manual_label_tmp$sample),function(x){
        tmp <- manual_label_tmp[manual_label_tmp$sample==x,]
        layer_region <- table(tmp$layer,tmp$region)
        layer_region <- as.data.frame(reshape2::melt(layer_region))
        layer_region$layer <- as.character(layer_region$Var1)
        layer_region$region <- as.character(layer_region$Var2)
        
        region_count <- table(tmp$region)
        layer_region$SR <- as.numeric(region_count[layer_region$region])
        layer_region$relative_value <- layer_region$value/layer_region$SR
        layer_region$Sample <- x
        layer_region <- layer_region[,c('Sample','layer','region','relative_value')]
        return(layer_region)}))
    
    Sample_plot_all <- Sample_plot_all[Sample_plot_all$region!='None',]
    IS <- c('GI','DI','AI')
    TE <- c('TE1','TE2','TE3')
    Mot <- c('A4ab','A4c','A6DC','A6DR','A6M')
    V1 <- c('V1')
    ACC <- c('A24a','A24b','A24c','A24d')
    Pir <- c('Pir')
    Ent <- c('Ent')
    OFC <- c('OPAl','OPro','ProM')
    
    Sample_plot_all$L2 <- 'Other'
    Sample_plot_all[Sample_plot_all$region%in%IS,'L2'] <- 'IS'
    Sample_plot_all[Sample_plot_all$region%in%Mot,'L2'] <- 'Mot'
    Sample_plot_all[Sample_plot_all$region%in%V1,'L2'] <- 'V1'
    Sample_plot_all[Sample_plot_all$region%in%TE,'L2'] <- 'TE'
    Sample_plot_all[,'L2'] <- factor(Sample_plot_all[,'L2'],levels = c('IS','TE','Other','Mot','V1'))
    
    Sample_plot_all$L3 <- 'Other'
    Sample_plot_all[Sample_plot_all$region%in%IS,'L3'] <- 'IS'
    Sample_plot_all[Sample_plot_all$region%in%Mot,'L3'] <- 'Mot'
    Sample_plot_all[Sample_plot_all$region%in%V1,'L3'] <- 'V1'
    Sample_plot_all[Sample_plot_all$region%in%TE,'L3'] <- 'TE'
    Sample_plot_all[,'L3'] <- factor(Sample_plot_all[,'L3'],levels = rev(c('IS','TE','Other','Mot','V1')))
    
    Sample_plot_all$L4 <- 'Other'
    Sample_plot_all[Sample_plot_all$region%in%ACC,'L4'] <- 'ACC'
    Sample_plot_all[Sample_plot_all$region%in%Mot,'L4'] <- 'Mot'
    Sample_plot_all[Sample_plot_all$region%in%V1,'L4'] <- 'V1'
    Sample_plot_all[,'L4'] <- factor(Sample_plot_all[,'L4'],levels = c('V1','Other','Mot','ACC'))
    
    Sample_plot_all$L2_Allo <- 'Other'
    Sample_plot_all[Sample_plot_all$region%in%Pir,'L2_Allo'] <- 'Pir'
    Sample_plot_all[Sample_plot_all$region%in%Ent,'L2_Allo'] <- 'Ent'
    Sample_plot_all[,'L2_Allo'] <- factor(Sample_plot_all[,'L2_Allo'],levels = c('Ent','Pir','Other'))

    my_comparisons1 <- list(c("Ent", "Other"), c("Pir", "Other"))
    my_comparisons2 <- list(c("IS", "Other"), c("TE", "Other"),  c("Mot", "Other"), c("V1", "Other"))
    my_comparisons3 <- list(c("IS", "Other"), c("TE", "Other"),  c("Mot", "Other"), c("V1", "Other"))
    my_comparisons4 <- list(c("V1", "Other"), c("Mot", "Other"), c("ACC", "Other"))
    my_comparisons5 <- list(c("V1", "Other"), c("Mot", "Other"), c("ACC", "Other"))
    my_comparisons6 <- list(c("OFC", "Other"), c("IS", "Other"), c("V1", "Other"))
    
    
    p1 <- ggplot(Sample_plot_all[Sample_plot_all$layer=='L2/3/4',],aes(x=L2_Allo,y=relative_value,color=layer))+
    geom_boxplot(fill=NA,lwd=1,outlier.shape = NA,show.legend = F,width=0.8)+ 
  stat_compare_means(comparisons = my_comparisons1)+
    scale_color_manual(values = layer_color,breaks = names(layer_color))+
    theme_classic()+
    theme(text = element_text(size=16))
    
    p2 <- ggplot(Sample_plot_all[Sample_plot_all$layer=='L2',],aes(x=L2,y=relative_value,color=layer))+
    geom_boxplot(fill=NA,lwd=1,outlier.shape = NA,show.legend = F,width=0.8)+
    stat_compare_means(comparisons = my_comparisons2)+
    scale_color_manual(values = layer_color,breaks = names(layer_color))+
    theme_classic()+
    theme(text = element_text(size=16))
    
    p3 <- ggplot(Sample_plot_all[Sample_plot_all$layer=='L3',],aes(x=L3,y=relative_value,color=layer))+
    geom_boxplot(fill=NA,lwd=1,outlier.shape = NA,show.legend = F,width=0.8)+
    stat_compare_means(comparisons = my_comparisons3)+
    scale_color_manual(values = layer_color,breaks = names(layer_color))+
    theme_classic()+
    theme(text = element_text(size=16))
    
    p4 <- ggplot(Sample_plot_all[Sample_plot_all$layer=='L4',],aes(x=L4,y=relative_value,color=layer))+
    geom_boxplot(fill=NA,lwd=1,outlier.shape = NA,show.legend = F,width=0.8)+
    stat_compare_means(comparisons = my_comparisons4)+
    scale_color_manual(values = layer_color,breaks = names(layer_color))+
    theme_classic()+
    theme(text = element_text(size=16))
    p <- p4+p1+p2+p3
    return(p)
}