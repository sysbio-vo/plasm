library(Seurat)
library(DoubletDecon)
library(dplyr)

# Load metadata

samples <- read.table("../meta/samples.tsv", sep="\t", stringsAsFactors = F , header = T)

#meta <- read.table("../meta/clusters_term_FH.tsv", sep="\t", stringsAsFactors = F, header = T)
#meta <- read.table("../meta/clusters_FH_villi.csv", sep=",", stringsAsFactors = F, header = T)
trims = c("first", "term")
tissues = c("villi", "placenta", "decidua")
trim = trims[2]
tissue = tissues[2]
diagnoses = c("PE", "C")
diagnosis = diagnoses[2]

# Load dataset
ids <- samples[which((samples$Tissue==tissue) & (samples$Trimester==trim) & (samples$Diagnosis==diagnosis)),]$ID
folders <- paste(ids, stringr::str_to_title(tissue), sep="-")

plc.data <- Read10X_h5(paste("../data/CR3/", folders[1], "/filtered_feature_bc_matrix.h5", sep=""))

plc <- CreateSeuratObject(counts = plc.data, min.cells = 3, min.features = 200, 
                           project = "plasm")

#cols <- colnames(plc@assays$RNA@counts)
#cols <- strsplit(cols, "-")
#cols <- sapply(cols, "[", 2)
#colnames(plc@meta.data)
#plc@meta.data$specimen <- cols
#meta <- meta[which(meta$Barcode %in% colnames(plc)),]

#plc@meta.data$clusters <- meta$Cluster[match(colnames(plc), meta$Barcode)]
#plc@meta.data$clusters[which(plc@meta.data$clusters=="")] <- "unknown"

### QC
plc[["percent.mt"]] <- PercentageFeatureSet(object = plc, pattern = "^MT-")
VlnPlot(object = plc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = plc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plc <- subset(x = plc, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 1)

### scTransform
plc <- SCTransform(object = plc, verbose = FALSE)
#plc <- SCTransform(object = plc, vars.to.regress = "percent.mt", verbose = FALSE)

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

DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "specimen", alpha=0.5)
        + NoLegend()
DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "specimen", alpha=0.5)
+ NoLegend()

DimPlot(object = plc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "clusters", alpha=0.5)
+ NoLegend()
DimPlot(object = plc, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "clusters_term_FH", alpha=0.5)
+ NoLegend()



FeaturePlot(object = plc, features = c("F13A1"))
FeaturePlot(object = plc, features = c("LYZ"))
FeaturePlot(object = plc, features = c("CD2"))
FeaturePlot(object = plc, features = c("KLRF1"))
FeaturePlot(object = plc, features = c("BACH2"))

FeaturePlot(object = plc, features = c("KLRF1"), reduction = "tsne")

FeaturePlot(object = plc, features = c("IGF2"))
FeaturePlot(object = plc, features = c("PAPPA2"))


head(plc@assays$RNA@counts)

data <- as.matrix(x = GetAssayData(plc, slot = "scale.data"))


data <- data[VariableFeatures(plc),]
data <- round(data, 2)
nrow(data)
ncol(data)
write.table(data, "scaled_counts.tsv", sep="\t", quote = F)
write.table(meta, "meta.tsv", sep="\t", quote = F, row.names = F)

genes <- rownames(data)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
ensg <- genes(EnsDb.Hsapiens.v86, filter=list(GeneNameFilter(genes),GeneIdFilter("ENSG", "startsWith")), 
      return.type="data.frame", columns=c("gene_id"))

ensg <- ensg[!duplicated(ensg$gene_name),]

data <- data[which(rownames(data) %in% ensg$gene_name), ]
nrow(data)
data <- data[match(ensg$gene_name, rownames(data)),]
head(rownames(data))
head(ensg)
rownames(data) <- ensg$gene_id

colnames(meta) <- c("Cell",	"cell_type")
meta <- meta[which(meta$cell_type!=""),]
meta$cell_type <- sapply(strsplit(meta$cell_type, " "), "[[", 2)
meta$cell_type[which(meta$cell_type=="SCT2")] <- "SCT"
meta$cell_type[which(meta$cell_type=="VCT2")] <- "VCT"
meta$cell_type[which(meta$cell_type=="EC1")] <- "EC"
meta$cell_type[which(meta$cell_type=="EC2")] <- "EC"
meta$cell_type[which(meta$cell_type=="HB1")] <- "HB"
meta$cell_type[which(meta$cell_type=="HB2")] <- "HB"


meta <- meta[which(meta$Cell %in% colnames(data)),]
data <- data[,which(colnames(data) %in% meta$Cell)]
ncol(data)
write.table(meta, "meta_villi.tsv", sep="\t", quote = F, row.names = F)
write.table(data, "scaled_counts_villi.tsv", sep="\t", quote = F)
