#Title: Cell Chat comparing 2 conditions 
#Author : Evgenia Ntaouka 
#Date: 10/5/2023

library(CellChat)
library(patchwork)
library(ComplexHeatmap)

setwd("C:/Users/zenia/Desktop/neurochemistry/CellChat/")

#Load the wildtype and mutant cellchat objects 

cellchat.MT <- readRDS("C:/Users/zenia/Desktop/neurochemistry/CellChat/CellChat_Full_Analysis_log_mutant.rds")
cellchat.WT <- readRDS("C:/Users/zenia/Desktop/neurochemistry/CellChat/CellChat_Full_Analysis_log_wildtype.rds")


#Subset only the cell groups you want 
cellchat.WT2 <- subsetCellChat(cellchat.WT, idents.use = c("OPCs", "Astrocytes", "Endothelial cells"))

cellchat.MT2<- subsetCellChat(cellchat.MT, idents.use = c("OPCs", "Astrocytes", "Endothelial cells"))

#Merge the cellchat objects
object.list <- list(WT = cellchat.WT2, MT = cellchat.MT2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


#Compare interactions 
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Subset communication 

Oligo_Astro_net_mutant <- subsetCommunication(cellchat)

Oligo_Astro_netP_mutant <- subsetCommunication(cellchat, slot.name = "netP") 
#To identify the interaction between which cell populations showing significant 
#changes, CellChat compares the number of interactions and interaction strength among different cell populations.
#Differential number of interactions or interaction strength among different cell populations


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#To do the heatmap with the right order of cells  use my netvisualheatmap2 

netVisual_heatmap2(cellchat.WT, color.heatmap = c("#F9CBD3","#7A3803"), measure = "weight", ) +  netVisual_heatmap2(cellchat.MT, color.heatmap = c("#F9CBD3","#7A3803"), measure = "weight")

#To visualize specific signalling pathways and the cell type interaction directly for mutant and wildtyep condition 

pathways.show <- c("NRXN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = c("#F9CBD3","#7A3803"),title.name = paste(pathways.show, "signaling ",names(object.list)[i]), measure = "weight")
}

#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))



#subset seurat to top 3 to plot genes

Idents(gaba_seurat) <- "level_4"
gaba_seurat2 <- subset(gaba_seurat, idents = c("Endothelial cells", "OPCs","Astrocytes"))

#Vlnplot from seurat package, Both for NRXN and WNT signallin gpathway
VlnPlot(gaba_seurat2,features = c("Nrxn1","Nrxn2","Nlgn3","Nlgn2","Nlgn1" ), group.by ="level_4", split.by = "mutation_status", pt.size = 0, cols =c( "#4b0082", "#ffee75"))


#Here i added some info in the function of netVisual_bubble so that i can have the exact numbers in the scale representing the communication 
#probability 
#Plot for the wildtype and mutant cellchat
net(cellchat.WT2, remove.isolate = FALSE, max.dataset = TRUE, signaling = c("WNT"), angle.x = 45, thresh = 0.05) +
net(cellchat.MT2, remove.isolate = FALSE, max.dataset = TRUE, signaling = c("WNT"), angle.x = 45, thresh = 0.05) 
