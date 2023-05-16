#Milo in fulldataset 22q11.2
# Zenia 7/11/22 



#Loading libraries
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

# define path
data_path <- "the paths for your file"
# load seurat object
data <- Connect(data_path) 
gaba_seurat <- as.Seurat(data, features = "var_names", cells = "obs_names")

#Normalize the dataset with SCtransform 
gaba_seurat <-SCTransform(gaba_seurat)

gaba_seurat <- RunPCA(gaba_seurat, npcs = 30)
gaba_seurat <- RunUMAP(gaba_seurat, dims = 1:30 , reduction = "pca", assay= SCT)

gaba_seurat<- RunTSNE(gaba_seurat)

#UMAP coloured by different slots from meta data
pdf("/results/Milo/UMAP_level4.pdf")
DimPlot(gaba_seurat, reduction = "umap", group.by = "level_4", label = TRUE)
dev.off()

pdf("results/Milo/UMAP_mutation_status.pdf")
DimPlot(gaba_seurat, reduction = "umap", group.by = "mutation_status")
dev.off()

pdf("/results/Milo/UMAP_sex.pdf")
DimPlot(gaba_seurat, reduction = "umap", group.by = "sex")
dev.off()

umap_plot_seurat <- UMAPPlot(gaba_seurat, group.by = "level_4", label=TRUE)


# transform  Seurat object to SingleCellExperiment object
gaba_sce <- as.SingleCellExperiment(gaba_seurat)

#Step 2: Visualize the data
umap_plot_seurat <- UMAPPlot(gaba_seurat, group.by = "level_4", label=TRUE)
umap_plot_seurat

#step 3 : Create Milo 


gaba_milo <- Milo(gaba_sce)
gaba_milo



#Step 4 : Build kNN Graph 
#k = S*3-5, rule of thumb , k it depends on the amount of cells
#k indeed, maybe something k~20 or so. With a data set that small you can play
#around a bit with the parameters too. For the proportion, I'd say ~0.3 is fine


gaba_milo <- buildGraph(gaba_milo, k =60, d = 30, reduced.dim = "PCA" ) #k>=45


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


gaba_milo <- makeNhoods(gaba_milo, prop = 0.05, k = 60 , d=30, refined = TRUE, reduced_dims = "PCA" )


pdf("/results/Milo/NhoodSizeHist.pdf")
plotNhoodSizeHist(gaba_milo)  
dev.off()

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

da_results_gaba <- testNhoods(gaba_milo, design = ~ mutation_status , design.df =gaba_design , reduced.dim = "PCA" ) 
head(da_results_gaba)

#corrected p values 

da_results_gaba %>%
  arrange(SpatialFDR) %>%
  head()

#Step 10 : Checking the resutls 

#We first inspect the distribution of uncorrected P values, to verify that the test was balanced.

pdf("/results/Milo/Pvalue.pdf")
ggplot(da_results_gaba, aes(PValue)) + geom_histogram(bins=50) #In R, the aes() function is often used within other graphing elements to specify the desired aesthetics. 
dev.off()

#Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).

pdf("/results/Milo/spatialfdr.pdf")
ggplot(da_results_gaba, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
dev.off()


#Build nhood graph 
gaba_milo <- buildNhoodGraph(gaba_milo)



#Plot Nhood grpah coloured by LogFC values 
nh_graph_pl_fdr <- plotNhoodGraphDA(gaba_milo, da_results_gaba,alpha=0.05, layout= "UMAP")

pdf("/results/Milo/umap_da_plot.pdf")
umap_plot_seurat + nh_graph_pl_fdr +
  plot_layout(guides="collect")
dev.off()


pdf("/results/Milo/MANeigh.pdf")
plotNhoodMA(da_results_gaba,alpha=0.05)
dev.off()

#To do this, we assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood. We can label neighbourhoods in the results data.frame using the function annotateNhoods.
#This also saves the fraction of cells harbouring the label.

da_results_gaba <- annotateNhoods(gaba_milo, da_results_gaba, coldata_col = "level_4")
head(da_results_gaba)

#While neighbourhoods tend to be homogeneous, we can define a threshold
#for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
pdf("/results/Milo/level4fraction.pdf")
ggplot(da_results_gaba, aes(level_4_fraction)) + geom_histogram(bins=50)
dev.off()


#Now we can visualize the distribution of DA Fold Changes in different cell types

pdf("/Milo/DAbeeswarm.pdf")
plotDAbeeswarm(da_results_gaba, group.by ="level_4", alpha = 0.05)
dev.off()




