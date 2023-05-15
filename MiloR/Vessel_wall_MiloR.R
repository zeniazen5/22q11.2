#DA with MiloR in Vessel wall cells 
# Zenia 7/11/22 



#Loading the libraries
library(SeuratDisk)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(MouseGastrulationData)
library(igraph)
library(RANN)




#If you have saved the seurat file then the part 1 is not needed 

## Part 1##

#The vessel wall cells are already normalized with SC transform 
setwd("C:/Users/zenia/Desktop/Code/Milo Zenia/Milo Vessel wall")
gaba_seurat<- LoadH5Seurat("Vessel wall SC .h5Seurat")

#Divide to have only endothelial cells 

Idents(gaba_seurat) <- "level_4"
gaba_seurat2<- subset(gaba_seurat, idents = "Endothelial cells")


#Umap of the endothelial cells 
umap_plot_seurat <- DimPlot(gaba_seurat3, reduction = "umap", group.by = "level_4", label = FALSE)+theme_void() +   scale_color_manual(breaks = cts, values = colors)   + NoAxes()


# transform gaba Seurat object to SingleCellExperiment object
gaba_sce <- as.SingleCellExperiment(gaba_seurat3)


#step 3 : Create Milo 


gaba_milo <- Milo(gaba_sce)
gaba_milo


#Step 4 : Build kNN Graph 
#k = S*3-5

gaba_milo <- buildGraph(gaba_milo, k =15, d = 30, reduced.dim = "PCA" , ) #k>=45

#Step 5 :  Define Neighborhoods

#For efficiency, we don’t test for DA in the neighbourhood of every cell, 
#but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015.
#
#prop: the proportion of cells to randomly sample to start with. We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster).
#refined: indicates whether you want to use the sampling refinement algorithm, or just pick cells at random. 
#The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal.

#we plot the distribution of neighbourhood sizes (i.e. how many cells form each neighbourhood) 
#to evaluate whether the value of k used for graph building was appropriate.

#Empirically, we found it’s best to have a distribution peaking above 20.
#Otherwise you might consider rerunning makeNhoods increasing k and/or prop.


gaba_milo <- makeNhoods(gaba_milo, prop = 0.2, k = 15 , d=30, refined = TRUE, reduced_dims = "PCA")

plotNhoodSizeHist(gaba_milo)  


#Step 6 : Calculate cells in neighborhoods


#how many cells frome each sample in each neighborhood 

gaba_milo <- countCells(gaba_milo, meta.data = data.frame(colData(gaba_milo)), sample="mouse_id")
head(nhoodCounts(gaba_milo))

#Step 7 : DA in neighborhoods , set experimental design

gaba_design <- data.frame(colData(gaba_milo))[,c("mouse_id", "mutation_status")]
gaba_design <- distinct(gaba_design) #used to select distinct or unique rows from the R data frame
rownames(gaba_design) <- gaba_design$mouse_id

gaba_design


#Step 8 : Calculate neighborhood connectivity 

#N.B. this step is the most time consuming of the analysis workflow and might 
#take a couple of minutes for large datasets).

gaba_milo <- calcNhoodDistance(gaba_milo, d=30, reduced.dim = "PCA")

#' This function calculates the mean expression of each feature in the
#' Milo object stored in the assays slot. Neighbourhood expression data are
#' stored in a new slot code nhoodExpression.
#'
gaba_milo <- calcNhoodExpression(gaba_milo)
dim(nhoodExpression(gaba_milo))


#Step 9 : Testing 


da_results_gaba <- testNhoods(gaba_milo, design = ~ mutation_status , design.df =gaba_design , reduced.dim = "PCA",fdr.weighting = "graph-overlap" ) 
head(da_results_gaba)



#corrected p values 

da_results_gaba %>%
  arrange(SpatialFDR) %>%
  head()

#Step 10 : Checking the resutls 

#We first inspect the distribution of uncorrected P values, to verify that the test was balanced.

ggplot(da_results_gaba, aes(PValue)) + geom_histogram(bins=50) #In R, the aes() function is often used within other graphing elements to specify the desired aesthetics. 


#Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).

ggplot(da_results_gaba, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

#Build Graph 
 
gaba_milo <- buildNhoodGraph(gaba_milo)


#LogFC DA graph 
nh_graph_pl_fdr <- plotNhoodGraphDA(gaba_milo, da_results_gaba,alpha=0.1, layout= "UMAP")



nh_graph_pl <- plotNhoodGraph(gaba_milo,   layout= "UMAP", )


umap_plot_seurat + nh_graph_pl_fdr +
  plot_layout(guides="collect")

#umap_plot_seurat + nh_graph_pl_Pvalue +
#plot_layout(guides="collect")


plotNhoodMA(da_results_gaba,alpha=0.05)


#To do this, we assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood. We can label neighbourhoods in the results data.frame using the function annotateNhoods.
#This also saves the fraction of cells harbouring the label.

da_results_gaba <- annotateNhoods(gaba_milo, da_results_gaba, coldata_col = "mutation_status")
head(da_results_gaba)

#While neighbourhoods tend to be homogeneous, we can define a threshold
#for celltype_fraction to exclude neighbourhoods that are a mix of cell types.


ggplot(da_results_gaba, aes(mutation_status_fraction)) + geom_histogram(bins=50)

#new :  to exclude neighbourhoods that are a mix of cell types. 

#da_results_gaba$mutation_status <- ifelse(da_results_gaba$mutation_status_fraction < 0.7, "Mixed", da_results_gaba$mutation_status)


#Now we can visualize the distribution of DA Fold Changes in different cell types

plotDAbeeswarm(da_results_gaba, group.by ="mutation_status", alpha = 0.1) +guides(shape=FALSE)




