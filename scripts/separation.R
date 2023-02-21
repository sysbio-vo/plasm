library(Seurat)
library(DoubletDecon)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load metadata
samples <- read.table("../meta/metadata.csv", sep=",", stringsAsFactors = F , header = T)
annotation <- read.table("../meta/Placenta_anno_semifinal_22012021.csv", sep=",", stringsAsFactors = F , header = T)

# Load dataset
i = 1
plc.data <- Read10X_h5(paste("../data/latest/", samples$h5name[i], sep=""))
plc <- CreateSeuratObject(counts = plc.data, min.cells = 3, min.features = 200, 
                          project = "plasm")

# QC
plc[["percent.mt"]] <- PercentageFeatureSet(object = plc, pattern = "^MT-")
VlnPlot(object = plc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plc <- subset(x = plc, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 1)

# Doublets removal
#newFiles=Improved_Seurat_Pre_Process(plc, num_genes=50, write_files=FALSE)

# Normalization
#plc <- SCTransform(object = plc, verbose = FALSE)

# Dimentionality reduction
plc <- RunPCA(object = plc, features = VariableFeatures(object = plc))

#plc <- JackStraw(object = plc, num.replicate = 100)
#plc <- ScoreJackStraw(object = plc, dims = 1:20)
#JackStrawPlot(object = plc, dims = 1:15)

# Clustering
plc <- FindNeighbors(object = plc, dims = 1:30)
plc <- FindClusters(object = plc, resolution = 0.5)

gc()
#plc <- RunTSNE(object = plc, dims.use = 1:10, do.fast = TRUE)
#DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
plc <- RunUMAP(object = plc, dims = 1:30)
DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Mito genes plot
plc@meta.data$mito <- 0
plc@meta.data[which(plc@meta.data$percent.mt>0.5),]$mito <- 1
DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "mito") + NoLegend()

# Sex genes plots

## Female
fplot <- FeaturePlot(object = plc, features = c("XIST"))

## Male
male_genes <- c("RPS4Y1", "EIF1AY1", "KDM5D", "DDX3Y", "DAZ", "BPY2", "SRY", "TSPY", "ZFY", "CDY", "TGIF2LY", "PRY", "VCY", "EIFAY")
features <- rownames(plc@assays$RNA@meta.features)
mgenes <- male_genes[male_genes %in% features]

# Calculate the difference between the average expression levels of a gene set and random control genes (all genes in our case)
# Positive numbers indicate higher expression level than random, i.e. it suggests that this module of genes is expressed in a particular cell
# more highly than would be expected, given the average expression of this module across the population.
plc <- AddModuleScore(object = plc, features = list(mgenes), name = "male_genes")
mplot <- FeaturePlot(object = plc, features = "male_genes1") + ggtitle(paste("Male Genes:", paste(mgenes, collapse = ", ")))
fplot + mplot

which(plc$male_genes1>0)
RidgePlot(plc, features = mgenes, ncol = 2)
RidgePlot(plc, features = "XIST", ncol = 2, group.by = "orig.ident")
BarPlot(plc, features = "XIST", ncol = 2, group.by = "orig.ident")
# Add cluster information
annotation$barcodes <- sapply(strsplit(annotation$X, "_"), "[", 2)
meta <- annotation[which(annotation$stim==samples$stim[i]),]
meta <- meta[which(meta$barcodes %in% colnames(plc)),]

plc@meta.data$clusters <- meta$cell_type_semifina[match(colnames(plc), meta$barcodes)]
plc@meta.data$clusters[which(plc@meta.data$clusters=="")] <- "unknown"

ctypes <- DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "clusters")

fplot + mplot + ctypes
# New clusters
pident=as.factor(meta$cell_type[match(colnames(plc), meta$barcodes)])
names(pident)=colnames(plc)  
plc@meta.data$orig.ident=pident


which(rownames(plc@assays$RNA@counts)=="XIST")

mean(plc@assays$RNA@counts[9697,])

max(plc@assays$RNA@counts)