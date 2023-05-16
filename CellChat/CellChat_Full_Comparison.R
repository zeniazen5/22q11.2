#CellChat analysis is LogTransform data Wt vs MT 
#7/3/2023

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

#Cell Chat libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)



 
#Load seurat 
setwd("/results/SCtransform")
gaba_no<- LoadH5Seurat("Seurat_Full.h5seurat")


#LogNormalizeData 
gaba_seurat <- NormalizeData(gaba_no, scale.factor = 10000 , normalization.method = "LogNormalize")


#Save the Log Normalize seurat 
setwd("/CellChat")
SaveH5Seurat(gaba_seurat, "Full_Log.h5Seurat")

#Make mean variance plot to check how your data look normalized
dg<- gaba_seurat@assays$RNA@data

gene_attr <- data.frame(mean = rowMeans(dg), var = apply(dg, 1, var))
gene_attr <- replace(gene_attr, gene_attr==0,1)
gene_attr$log_mean <- log10(gene_attr$mean)


#Save the plot 

png("results/CellChat/mean_variance.png")
ggplot( gene_attr, aes( log_mean,var )) + geom_point() +labs(x= "Log-mean of log_normalized data", y= "Variance of log_normalized data", title = "Mean-Variance plot") +
  theme_classic()
dev.off()


#Divide wildtype to do the separate cellchat analysis 
Idents(gaba_seurat)<-"mutation_status"
gaba_seurat2 <- subset(gaba_seurat, idents = c("wildtype"))

#Double check if you have only the mutant status
unique(gaba_seurat2@meta.data[["mutation_status"]])

#Save the wildtype seurat file
setwd("/results/CellChat")
SaveH5Seurat(gaba_seurat2, "Full_Log_wildtype.h5Seurat")

#CellChat 

########----------- CELL CHAT ANALYSIS FOR WILDTYPE--------#########
#Create cellChat object 

gaba_meta <- gaba_seurat2@meta.data
cellchat <- createCellChat(object = gaba_seurat2, meta = gaba_meta , group.by = "level_4")

# If meta slot is empty 
#cellchat <- addMeta(cellchat, meta = meta)
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#levels(cellchat@idents) # show factor levels of the cell labels
#groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


#Add the database for the cellchat analysis

CellChatDB <- CellChatDB.mouse # use CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB, you can also use a subset 

# set the used database in the object
cellchat@DB <- CellChatDB.use

#Save the cell chat object

setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_beforeOEG_log_wildtype.rds")

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

#Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Save the object after identification of overexpressed genes
setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_afterOEG_log_wildtype.rds")

#Compute the Communication probability: the "triMean" was used meaning that  
#implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%.
#The population size = TRUE to take into account the different size of cell types, do.fast= TRUE and nboot =20 in order to reduce the memory capacity
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean", population.size = TRUE, do.fast = TRUE, nboot = 20 ) #do.fast = true, nboo = less permutation
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
Net <- subsetCommunication(cellchat) 

NetP<- subsetCommunication(cellchat, slot.name = "netP") 
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


#CellChat computes the communication probability on signaling pathway level by summarizing 
#the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.

cellchat <- computeCommunProbPathway(cellchat)

#Save Cell chat 

setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_Analysis_log_wildtype.rds")

#Export dataframe in Excel with all the LR pairs and signalling pathways 
library("writexl")
write_xlsx(Net,"/results/CellChat/Net_Full_wildtype.xlsx")

write_xlsx(NetP,"/results/CellChat/NetP_Full_wildtype.xlsx")

########----------- CELL CHAT ANALYSIS FOR MUTANT--------#########
rm(cellchat)
rm(gaba_meta)
#Divide wildtype to do the separate cellchat analysis 
Idents(gaba_seurat)<-"mutation_status"
gaba_seurat3 <- subset(gaba_seurat, idents = c("mutant"))

unique(gaba_seurat3@meta.data[["mutation_status"]])

setwd("/results/CellChat")
SaveH5Seurat(gaba_seurat3, "Full_Log_mutant.h5Seurat")

#CellChat 


#Create cellChat object 

gaba_meta <- gaba_seurat3@meta.data
cellchat <- createCellChat(object = gaba_seurat3, meta = gaba_meta , group.by = "level_4")

# If meta slot is empty 
#cellchat <- addMeta(cellchat, meta = meta)
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#levels(cellchat@idents) # show factor levels of the cell labels
#groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


CellChatDB <- CellChatDB.mouse # use CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB, you can also use a subset 

# set the used database in the object
cellchat@DB <- CellChatDB.use


setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_beforeOEG_log_mutant.rds")

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_afterOEG_log_mutant.rds")

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean", population.size = TRUE, do.fast = TRUE, nboot = 20 ) #do.fast = true, nboo = less permutation
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
Net <- subsetCommunication(cellchat) 

NetP<- subsetCommunication(cellchat, slot.name = "netP") 
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#CellChat computes the communication probability on signaling pathway level by summarizing 
#the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.

cellchat <- computeCommunProbPathway(cellchat)

#Save Cell chat 

setwd("/results/CellChat")

saveRDS(cellchat, file = "CellChat_Full_Analysis_log_mutant.rds")

#Export dataframe in Excel 
library("writexl")
write_xlsx(Net,"/results/CellChat/Net_Full_mutant.xlsx")

write_xlsx(NetP,"/results/CellChat/NetP_Full_mutant.xlsx")





