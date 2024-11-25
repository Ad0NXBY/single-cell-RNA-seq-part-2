#Install and load the required R libraries for this analysis
install.packages("dplyr")
install.packages("Seurat")
install.packages("patchwork")
install.packages("ggplot2")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Load in the .rds objects into the R environment
pbmc <- readRDS("pbmc3k_QCB_Session5(1).rds")

#Creating additional QC metrics to visualize==============================================
#Amount of Ribosomal Genes
#Calculating the percentage ribosomal content of the RNA for each barcode and variable put into pbmc@meta.data.slot
pbmc@meta.data$pt.Ribosomal <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
#Percentage of Largest gene
#Get the count for the largest gene per cell- The index position of the gene with the largest count- The name of the most highly expressed gene per cell
#Calculates the largest gene and its percentaged will be stashed into pbmc@meta.data
pbmc.nomalat <- pbmc[rownames(pbmc) != "MALAT1", ] #grab this out as its a non-coding nuclear gene expressed at high level, has big effect so measure separately
pbmc.nomalat$largest_count <- apply(pbmc.nomalat@assays$RNA@counts, 2, max)
pbmc.nomalat$largest_index <- apply(pbmc.nomalat@assays$RNA@counts, 2, which.max)
pbmc.nomalat$largest_gene <- rownames(pbmc.nomalat)[pbmc.nomalat$largest_index]
pbmc.nomalat$percent.Largest.Gene <- 100 * pbmc.nomalat$largest_count/pbmc.nomalat$nCount_RNA
pbmc@meta.data$largest_gene <- pbmc.nomalat$largest_gene
pbmc@meta.data$percent.Largest.Gene <- pbmc.nomalat$percent.Largest.Gene

#Clean up and delete pbmc.nomalat
rm(pbmc.nomalat)

#Visualize the new QC metrics 
VlnPlot(pbmc, features = c("pt.Ribosomal", "percent.Largest.Gene"))

#Explotation of the largest gene and its percentages across cell barcodes
I.genes <- as.data.frame(cbind(as.numeric(pbmc@meta.data$percent.Largest.Gene), pbmc@meta.data$largest_gene))
I.genes.sorted <- I.genes[order(I.genes$V1, decreasing = T), ]

#Perform linear dimensional reduction=======================================================
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
str(pbmc)

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Examine and visualize PCA results a few different ways
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#The plot is showing the top genes associated with reduction components

# Examine and visualize PCA results a few different ways. DimPlot produces a plot of the
# cell.embeddings for the first two PCs (PC1 and PC2).
DimPlot(pbmc, reduction = "pca")

# We can use the group.by option to colour by any other metadata column, today will will use
# largest_gene. We can also add labels to the plot. Finally we can add a call to the
# NoLegend() function to suppress the automatic colour legend which is drawn.

DimPlot(pbmc, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 3, repel = T) +
  NoLegend()

DimPlot(pbmc, reduction = "pca", label = TRUE, label.size = 3, repel = T) + NoLegend()

# Now examine the next two PCs (PC_3 and PC_4)
DimPlot(pbmc, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 2, repel = T,
        dims = 3:4) + NoLegend()
#Try 5:6, what do you see? It becomes a cluster with no distinct subgroups, making it difficult to distinguish between groups

#PCA plots limitations=========================================================================
#We can see that not all the useful information is captured int he first two prinicipal components, so how far down the set of PCs do we need to go to 
#Capture all the biologically relevant information, prior to clustering and visualization
#tiered strategy that uses heuristic, visualization and statiscal rresampling approaches to select most informative PCs

#1st method - Examine variance explained by the PCA analysis, by ranking stdev values for each PC calculated by RunPCA

# Extract the variance (stdev column) captured by each PC (PCs column), examine and plot the
# data
PCA_stdev <- as.data.frame(cbind(colnames(pbmc@reductions$pca@feature.loadings), pbmc@reductions$pca@stdev))
colnames(PCA_stdev) <- c("PCs", "stdev")

plot(rownames(PCA_stdev), PCA_stdev$stdev) #Creates a scree plot

#Use it using Seurat by using ElbowPlot function
x <- 1/mean(1/as.numeric(PCA_stdev$stdev))  # Calculates the harmonic mean of stdev for the first 50 PCs. We are using this to highlight the 'elbow' in the plot, what other value could you use?
ElbowPlot(pbmc, ndims = 50) + geom_hline(yintercept = x, color = "grey")

#2nd method - Visualization using Feature Heatmaps ordered according to their PCA scores.
#In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. 
#Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. 
#Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.

# Explore the first PC
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Explore the first 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#3rd method - perform a resampling test (The Jackstraw procedure) 
# NOTE: This process can take a long time for bigger datasets. More approximate techniques
# such as those implemented in ElbowPlot() can be used to reduce computation time.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# Plot the results of the JackStraw analysis
JackStrawPlot(pbmc, dims = 1:15)
#PC1-13 are significant due to their P-value being less than 0.05