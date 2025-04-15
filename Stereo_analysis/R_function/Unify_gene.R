#Homologous gene table from TableS2E
species_gene <- read.csv('Tabledata/mart_export.humanMacaqeMarmosetMouse.oneToOneOrth.ensembl91.20220428.txt',sep='\t')
change_gene <- function(genes,species1,species2){
    tmp_df <- species_gene
    rownames(tmp_df) <- tmp_df[,paste0(species1,'Gene')]
    return(tmp_df[genes,paste0(species2,'Gene')])}

#Search for genes of different species within the homolog list and standardize them to marmoset genes.
Find_homogene <- function(obj,species){
    obj[['RNA']]@meta.features$homo <- rownames(obj[['RNA']])%in%species_gene[,paste0(species,'Gene')]
    obj[['RNA']]@meta.features$homoname <- 'None'
    obj[['RNA']]@meta.features[obj[['RNA']]@meta.features$homo,]$homoname <- 
    change_gene(rownames(obj)[obj[['RNA']]@meta.features$homo],species,'marmoset')
    return(obj)
}
#Find the intersection gene in input list
inter_multi <- function(vector_list){
    inter <- vector_list[[1]]
    for(i in 2:length(vector_list)){
        inter <- intersect(inter,vector_list[[i]])
    }
    return(inter)}