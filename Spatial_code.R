#!/usr/bin/env Rscript

#Based on Vignette: https://satijalab.org/seurat/articles/spatial_vignette.html

library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(purrr)
library(topGO)
library(rgraphviz)
library(SingleCellExperiment)
library(BayesSpace)
library(scran)
library(scater)

# Load Data
tissue <- Load10X_Spatial(data.dir = "~/Desktop/Comparative_10X/Liver/FemaleB/outs/",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "B1_Tissue",
                          filter.matrix = TRUE,
                          to.upper = FALSE,
                          image = NULL)

# Filtering
selected_c <- WhichCells(tissue, expression = nFeature_Spatial > 100)
selected_f <- rownames(tissue)[Matrix::rowSums(tissue) > 3]

tissue <- subset(tissue, features = selected_f, cells = selected_c)

# Transform
tissue <- SCTransform(tissue, assay = 'Spatial', verbose = FALSE)

p2 <- VlnPlot(tissue, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- p2 + xlab( "Male A Spots") + ggtitle("UMIs Per Spot") +ylab("UMIs") + theme(axis.text.x=element_blank(),
                                                                                  axis.ticks.x=element_blank())

p3 <- VlnPlot(tissue, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p3 <- p3 + xlab( "Male A Spots") + ggtitle("Genes Per Spot") +ylab("Genes") + theme(axis.text.x=element_blank(),
                                                                                    axis.ticks.x=element_blank())
# Clustering
tissue <- RunPCA(tissue, assay = "SCT", verbose = FALSE)
tissue <- FindNeighbors(tissue, reduction = "pca", dims = 1:20)
tissue <- FindClusters(tissue, graph.name = "SCT_snn", verbose = FALSE)
tissue <- RunUMAP(tissue, reduction = "pca", dims = 1:20)
tissue <- RunTSNE(tissue, reduction = "pca", dims = 1:20)



new.cluster.ids <- c('Portal','Portal','Midzonal','Central','Erythro','Central','Non-Parenchymal')
names(new.cluster.ids) <- levels(tissue)
tissue <- RenameIdents(tissue, new.cluster.ids)

p1 <- SpatialDimPlot(tissue, label = TRUE, label.size = 3)


# Check Clusters
plotPCA <- DimPlot(tissue, reduction = "pca", label = TRUE, label.size = 5)
plotPCA
plotUMAP <- DimPlot(tissue, reduction = "umap", label = TRUE)
plotUMAP
plotTSNE <- DimPlot(tissue, reduction = 'tsne', label = TRUE)
plotTSNE


# Find Variable Features
tissue <- FindSpatiallyVariableFeatures(tissue, assay = "SCT", features = VariableFeatures(tissue)[1:1000], selection.method = "moransi")
morans_features <- SpatiallyVariableFeatures(tissue, selection.method = "moransi")

top_features <- head(morans_features, 9)
plotMorans <- SpatialFeaturePlot(tissue, features = top_features)
plotMorans2 <- VlnPlot(tissue, features = top_features) + theme(strip.text.x = element_blank())
plotMorans3 <- FeaturePlot(tissue, features = top_features)


# Plotting 

p4 <- SpatialFeaturePlot(tissue, features = "nCount_Spatial") + theme(legend.position = "right")
p5 <- SpatialFeaturePlot(tissue, features = "nFeature_Spatial") + theme(legend.position = "right")


#Spatial UMAP plot with each cluster separated, choose cluster with idents
p6 <- SpatialDimPlot(tissue, cells.highlight = CellsByIdentities(object = tissue),
                     facet.highlight = TRUE,
                     cols.highlight = c("#E87D72", "grey50"),
                     ncol = 4)


#Violin plot shows expression levels of top spatially variable features
p7 <- VlnPlot(tissue, features = top_features) + theme(strip.text.x = element_blank())

p8 <-  SpatialFeaturePlot(tissue, features = top_features, ncol = 3)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
markers <- FindAllMarkers(tissue, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Number of differentially expressed genes per cluster
PortalCount <- sum(markers$cluster == 'Portal')
MidzonalCount <- sum(markers$cluster == 'Midzonal')
CentralCount <- sum(markers$cluster == 'Central')
ErythroCount <- sum(markers$cluster == 'Erythro')
ImmuneCount <- sum(markers$cluster == 'Non-Parenchymal')

HeatMap <- DoHeatmap(tissue, features = top10$gene)

# Canonical Marker Genes for Cell Type Identity
marker_genes1 <- c('Acly', 'LOC102924223', 'LOC102928708', 'Glul','Jchain','LOC102925039') #From Literature
marker_genes2 <- c('Aldh1b1', 'EI127-mgr01', 'Glul', 'LOC102929018', 'LOC102913793') #Top per cluster
marker_genes3 <- c("LOC102929018", "LOC102928917", "LOC102928708", "LOC102902798") #Erythro
marker_genes4 <- c('Fxyd2','Col3a1','Acta2','LOC102913793','LOC102919474', 'Jchain') #Non Parenchymal
marker_genes5 <- c('Wnt2','Wnt9b','Rspo3','Dkk3') #Endothelial cell markers
marker_genes6 <- c('LOC102925039', 'Hal', 'Acly', 'Aldh1b1', 'Hsd17b13') # Portal
marker_genes7 <- c('Glul','LOC102924473', 'Oat', 'Lect2','Slc1a2') # Central
marker_genes8 <- c('Hint1', 'Cox7c', 'Apoc1', 'Fabp1', 'Mt2a', 'Mt1g', 'Ndufb1') #Midzonal
marker_genes9 <- c('LOC102924223', 'LOC102907986','LOC102908284','LOC102924473','LOC102909592','Tst') #xenobiotic cytochrome and others
marker_genes10 <- c('Hilpda','Plin2','Upp2','Serpina6') # Lipid genes of interest marker
marker_genes11 <- c('Itm2b','Apoa4','EI127-mgr01','LOC102908284', 'LOC102929018','Glul','LOC102913793') #Pre-annotated top DE gene in 7 clusters
marker_genes12 <- c('Tst','LOC102907986','LOC102906854') # more xenobiotic genes
marker_genes13 <- c('Acaca','Acacb', 'Elovl6') # lipogenesis genes
marker_genes14 <- c('Slc2a3', 'Slc2a1', 'Ldha') # Glycolysis genes

markervln <- VlnPlot(tissue, features = marker_genes1) + theme(strip.text = element_blank())


# TopGO analysis

# Cluster0
cluster <- subset(tissue, idents = 'Portal')

expr <- as.matrix(GetAssayData(cluster))
expressed.genesP <- markers$gene[markers$cluster == 'Portal']
all.genes <- rownames(expr)
writeLines(expressed.genesP, "clusterP_genes.txt")


# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genesP, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # Biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol") 

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Saves a new Go table after each iteration
df <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
filename <- 'ClusterP_GoTable.tsv'
write.table(df, file=filename, sep = '\t', row.names = F)


# Cluster Midzonal
cluster <- subset(tissue, idents = 'Midzonal')

expr <- as.matrix(GetAssayData(cluster))

expressed.genesM <- markers$gene[markers$cluster == 'Midzonal']
all.genes <- rownames(expr)
writeLines(expressed.genesM, "clusterM_genes.txt")


# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genesM, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # Biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol") 

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Saves a new Go table after each iteration
df <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
filename <- 'ClusterM_GoTable.tsv'
write.table(df, file=filename, sep = '\t', row.names = F)

# Cluster2
cluster <- subset(tissue, idents = 'Central')

expr <- as.matrix(GetAssayData(cluster))

expressed.genesC <- markers$gene[markers$cluster == 'Central']
all.genes <- rownames(expr)
writeLines(expressed.genesC, "clusterC_genes.txt")


# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genesC, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # Biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol") 

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Saves a new Go table after each iteration
df <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
filename <- 'ClusterC_GoTable.tsv'
write.table(df, file=filename, sep = '\t', row.names = F)


# Cluster Erythro
cluster <- subset(tissue, idents = 'Erythro')

expr <- as.matrix(GetAssayData(cluster))

expressed.genesE <- markers$gene[markers$cluster == 'Erythro']
all.genes <- rownames(expr)
writeLines(expressed.genesE, "clusterE_genes.txt")


# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genesE, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # Biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol") 

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Saves a new Go table after each iteration
df <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
filename <- 'ClusterE_GoTable.tsv'
write.table(df, file=filename, sep = '\t', row.names = F)


# Cluster Non-Parenchymal
cluster <- subset(tissue, idents = 'Non-Parenchymal')

expr <- as.matrix(GetAssayData(cluster))

expressed.genesN <- markers$gene[markers$cluster == 'Non-Parenchymal']
all.genes <- rownames(expr)
writeLines(expressed.genesN, "clusterN_genes.txt")


# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genesN, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # Biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol") 

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Saves a new Go table after each iteration
df <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
filename <- 'ClusterN_GoTable.tsv'
write.table(df, file=filename, sep = '\t', row.names = F)


# Number of spots per cluster
n_spots <- FetchData(tissue, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# overall number of features

cluster.id <- c("Portal","Midzonal",'Central',"Eryhtro", "Non-Parenchymal")
cluster.set <- unique(tissue$SCT_snn_res.0.8) # Factor with 7 levels
tissue$cluster_ids <- new.cluster.ids

tissue.1 <- subset(tissue, idents='Midzonal')

cells.c <- WhichCells(object = tissue, expression = SCT_snn_res.0.8 == 6)
sum(rowSums(tissue[['Spatial']]@counts[, cells.c ]) != 0)

# number features for each cluster
sapply(X = cluster.set, function(c) {
  cells.c <- WhichCells(object = tissue, expression = SCT_snn_res.0.8 == c)
  nFeature.c <- sum(rowSums(tissue[['Spatial']]@counts[, cells.c ]) != 0)
  return(nFeature.c)
})




### Load Data for BayesSpace ###

tissue <- Load10X_Spatial(data.dir = "~/Desktop/Comparative_10X/Liver/FemaleB/outs/",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "B1_tissue",
                          filter.matrix = TRUE,
                          to.upper = FALSE,
                          image = NULL)

# Filter
selected_c <- WhichCells(tissue, expression = nFeature_Spatial > 100)
selected_f <- rownames(tissue)[Matrix::rowSums(tissue) > 3]

tissue <- subset(tissue, features = selected_f, cells = selected_c)

# Normalization
tissue <- SCTransform(tissue, assay = 'Spatial', verbose = FALSE)
tissue <- DietSeurat(tissue, graphs = "pca")

# Convert to SCE
bs_tissue <- as.SingleCellExperiment(tissue)
colData(bs_tissue) = cbind(colData(bs_tissue), tissue@images$B1_tissue@coordinates)

bs_tissue <- spatialPreprocess(bs_tissue, platform="Visium",
                               n.PCs = 20, n.HVGs=2000, log.normalize=FALSE)

# Clustering

bs_tissue <- spatialCluster(bs_tissue, q=5, platform="Visium", d=20,
                            init.method="mclust", model="t", gamma=2,
                            nrep=20000, burn.in=100,
                            save.chain=TRUE)

plotClusters <- clusterPlot(bs_tissue)

bs_tissue.enhanced <- spatialEnhance(bs_tissue, q=5, platform="Visium", d=20,
                                     model="t", gamma=3,
                                     jitter_prior=0.3, jitter_scale=3.5,
                                     nrep=1000, burn.in=100,
                                     save.chain=TRUE)

plotEnhancedClusters <- clusterPlot(bs_tissue.enhanced)

# Cluster Marker Genes
markers <- list()
markers[["Portal Hepatocytes"]] <- c('LOC102925039', 'Hal', 'Acly')
markers[["Midzonal Hepatocyte"]] <- c('Enho','Clec4f')
markers[["Central Hepatocyte"]] <- c('Glul','LOC102924473', 'Oat', 'Lect2')
markers[["Erythro"]] <- c("LOC102929018", "LOC102928917", "LOC102928708", "LOC102902798")
markers[["Non-parenchymal"]] <- c('Col3a1','Acta2','Vim')

# Lipid metabolism genes
goi <- list()
goi[["Hilpda"]] <- c('Hilpda')
goi[['Plin2']] <- c('Plin2')
goi[['Upp2']] <- c('Upp2')
goi[['Serpina6']] <- c('Serpina6')

# Xenobiotic genes
goi2 <- list()
goi2[['CYP2E1']] <- c('LOC102924223')
goi2[['CYP1A1']] <- c('LOC102907986')
goi2[['CYP1A2']] <- c('LOC102908284')
goi2[['CYP2A3']] <- c('LOC102924473')
goi2[['CYP2B1']] <- c('LOC102909592')
goi2[['Rhodanese']] <- c('Tst')

bs_tissue.enhanced <- enhanceFeatures(bs_tissue.enhanced, bs_tissue,
                                      model="xgboost",
                                      feature_names=purrr::reduce(goi2, c),
                                      nrounds=0)

sum_counts <- function(bs_tissue, features) {
  if (length(features) > 1) {
    colSums(logcounts(bs_tissue)[features, ])
  } else {
    logcounts(bs_tissue)[features, ]
  }
}

# spot_expr <- purrr::map(markers, function(x) sum_counts(bs_tissue, x))
enhanced_expr <- purrr::map(goi2, function(x) sum_counts(bs_tissue.enhanced, x))

plot_expression <- function(bs_tissue, expr, name) {
  featurePlot(bs_tissue, expr, color=NA) +
    viridis::scale_fill_viridis(option="A") +
    labs(title=name, fill="Log-normalized\nexpression")
}


enhanced_plots <- purrr::imap(enhanced_expr, function(x, y) plot_expression(bs_tissue.enhanced, x, y))
patchwork::wrap_plots(enhanced_plots, ncol=3)
