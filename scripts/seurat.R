library(Seurat)
library(dplyr)

# Load the PBMC dataset
plc.data <- Read10X(data.dir = "../data/cellranger_short/Aggr-Term/")
#plc.data <- Read10X(data.dir = "../data/cellranger_short/Aggr-Villi/")
#plc.data <- Read10X(data.dir = "../data/cellranger_short/pbmc/")

plc <- CreateSeuratObject(counts = plc.data, min.cells = 3, min.features = 200, 
                           project = "plasm")


### QC

plc[["percent.mt"]] <- PercentageFeatureSet(object = plc, pattern = "^MT-")
VlnPlot(object = plc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plc <- subset(x = plc, subset = nFeature_RNA > 200 & nFeature_RNA < 2000)
### scTransform
#plc <- SCTransform(object = plc, verbose = FALSE)
plc <- SCTransform(object = plc, vars.to.regress = "percent.mt", verbose = FALSE)
### Filtering
#plc <- subset(x = plc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### Normalization
#plc <- NormalizeData(object = plc, normalization.method = "LogNormalize", scale.factor = 10000)

### Variable genes
#plc <- FindVariableFeatures(object = plc, selection.method = "vst", nfeatures = 2000)

### Scale
#all.genes <- rownames(x = plc)
#plc <- ScaleData(object = plc, features = all.genes)

### Dimentionality reduction
plc <- RunPCA(object = plc, features = VariableFeatures(object = plc))

plc <- JackStraw(object = plc, num.replicate = 100)
plc <- ScoreJackStraw(object = plc, dims = 1:20)
JackStrawPlot(object = plc, dims = 1:15)

### Clustering
plc <- FindNeighbors(object = plc, dims = 1:30)
plc <- FindClusters(object = plc, resolution = 0.5)

gc()
plc <- RunTSNE(object = plc, dims.use = 1:10, do.fast = TRUE)
plc <- RunUMAP(object = plc, dims = 1:30)
DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

##
FeaturePlot(object = plc.tsne, features.plot = c("LYZ"), cols.use = c("grey", "blue"), 
                              reduction.use = "tsne")
FeaturePlot(object = plc.tsne, features.plot = c("THEMIS"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

FeaturePlot(object = plc.tsne, features.plot = c("SPTA1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
FeaturePlot(object = plc.tsne, features.plot = c("DPYD"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
##

plc@meta.data$mito <- 0
plc@meta.data[which(plc@meta.data$percent.mt>1),]$mito <- 1

DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "mito", alpha=0.5)
        + NoLegend()
DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "mito", alpha=0.5)
        + NoLegend()


plc@meta.data$doubl <- 0
plc@meta.data[which(plc@meta.data$nFeature_RNA>1000),]$doubl <- 1

DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "doubl", alpha=0.5)
+ NoLegend()
DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "doubl", alpha=0.5)
+ NoLegend()
