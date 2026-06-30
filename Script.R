#!/usr/bin/env Rscript

# Definir limite estrito de memória para o R (ajuste para os seus ~60GB disponíveis)
Sys.setenv(R_MAX_VSIZE = "60Gb")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(SingleR)
  library(celldex)
  library(SummarizedExperiment)
  library(pheatmap)
})

# =========================
# Paths
# =========================
base_dir <- "." # Raiz do repositório clonado
barcodes_file <- file.path(base_dir, "data", "GSE236581_barcodes.tsv.gz")
features_file <- file.path(base_dir, "data", "GSE236581_features.tsv.gz")
mtx_file      <- file.path(base_dir, "data", "GSE236581_counts.mtx.gz")
metadata_file <- file.path(base_dir, "data", "GSE236581_CRC-ICB_metadata.txt.gz")
outdir        <- file.path(base_dir, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Read matrix
# =========================
message("Reading matrix...")
mat <- Matrix::readMM(mtx_file)
barcodes <- fread(barcodes_file, header = FALSE)
features <- fread(features_file, header = FALSE)
meta <- fread(metadata_file, data.table = FALSE)

if (ncol(features) >= 2) {
  gene_names <- make.unique(as.character(features[[2]]))
} else {
  gene_names <- make.unique(as.character(features[[1]]))
}
cell_names <- as.character(barcodes[[1]])
rownames(mat) <- gene_names
colnames(mat) <- cell_names

# =========================
# Build Seurat object
# =========================
message("Creating Seurat object...")
obj <- CreateSeuratObject(
  counts = mat,
  project = "GSE236581",
  min.cells = 5,       
  min.features = 300    
)

rm(mat); gc()

# =========================
# Metadata integration
# =========================
meta <- as.data.frame(meta)
candidate_barcode_cols <- c("barcode", "Barcode", "cell", "Cell", "cell_id", "CellID")
barcode_col <- candidate_barcode_cols[candidate_barcode_cols %in% colnames(meta)][1]
if (is.na(barcode_col)) {
  barcode_col <- colnames(meta)[1]
}
meta[[barcode_col]] <- as.character(meta[[barcode_col]])
rownames(meta) <- meta[[barcode_col]]
common_cells <- intersect(colnames(obj), rownames(meta))
obj <- subset(obj, cells = common_cells)
meta <- meta[colnames(obj), , drop = FALSE]
obj <- AddMetaData(obj, metadata = meta)

write.csv(obj@meta.data, file.path(outdir, "metadata.csv"), row.names = TRUE)

# =========================
# QC metrics & Filtering
# =========================
message("Calculating QC metrics...")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

# Gerar gráficos de QC em PNG antes do filtro
png(file.path(outdir, "QC_violin.png"), width = 12, height = 5, units = "in", res = 300)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
dev.off()

png(file.path(outdir, "QC_VlnPlots.png"), width = 12, height = 5, units = "in", res = 300)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0.1))
dev.off()

png(file.path(outdir, "QC_Scatter.png"), width = 10, height = 5, units = "in", res = 300)
print(
  FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
)
dev.off()

message("Filtering cells...")
obj <- subset(
  obj,
  subset = nFeature_RNA >= 400 & 
           nFeature_RNA <= 6000 & 
           nCount_RNA <= 40000 & 
           percent.mt <= 15
)
gc()

# =========================
# DOWNSAMPLING CRÍTICO (Para rodar em 64GB)
# =========================
message("Dataset muito grande para 64GB. Aplicando downsampling para 50.000 células...")
set.seed(42) 
if (ncol(obj) > 50000) {
  celulas_selecionadas <- sample(colnames(obj), size = 50000, replace = FALSE)
  obj <- subset(obj, cells = celulas_selecionadas)
}
gc() 

# =========================
# Standard Seurat workflow
# =========================
message("Running Standard Seurat Workflow...")

obj <- NormalizeData(obj)
gc()

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
gc()

message("Scaling variable features...")
features_to_scale <- VariableFeatures(obj)

obj <- ScaleData(
  obj, 
  features = features_to_scale,
  do.center = TRUE,
  do.scale = TRUE
)
gc()

obj <- RunPCA(obj, npcs = 50, features = features_to_scale)
gc()

obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.6)
obj <- RunUMAP(obj, dims = 1:30)

saveRDS(obj, file.path(outdir, "Seurat_GSE236581.rds"))

# =========================
# UMAP plots (PNG)
# =========================
png(file.path(outdir, "UMAP_clusters.png"), width = 8, height = 6, units = "in", res = 300)
print(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

png(file.path(outdir, "UMAP_Final.png"), width = 8, height = 6, units = "in", res = 300)
print(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

meta_cols <- colnames(obj@meta.data)
candidate_blueprint <- c("sample_type", "SampleType", "tissue", "Tissue", "site", "Site", "source", "Source")
blueprint_col <- candidate_blueprint[candidate_blueprint %in% meta_cols][1]
if (!is.na(blueprint_col)) {
  png(file.path(outdir, "UMAP_Blueprint.png"), width = 8, height = 6, units = "in", res = 300)
  print(DimPlot(obj, reduction = "umap", group.by = blueprint_col))
  dev.off()
}

umap_df <- Embeddings(obj, "umap") %>% as.data.frame()
umap_df$cell <- rownames(umap_df)
umap_df$cluster <- obj$seurat_clusters
write.csv(umap_df, file.path(outdir, "umap_coords.csv"), row.names = FALSE)

# =========================
# Cluster markers
# =========================
message("Finding cluster markers...")
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, file.path(outdir, "cluster_markers.csv"), row.names = FALSE)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

png(file.path(outdir, "Heatmap_Top5_Clusters.png"), width = 10, height = 12, units = "in", res = 300)
print(DoHeatmap(obj, features = unique(top5$gene), group.by = "seurat_clusters") + NoLegend())
dev.off()

avg_exp <- AverageExpression(obj, return.seurat = FALSE)$RNA
top_heat_genes <- unique(top5$gene)
avg_exp_sub <- avg_exp[intersect(top_heat_genes, rownames(avg_exp)), , drop = FALSE]

png(file.path(outdir, "Heatmap_Media_Genes.png"), width = 8, height = 10, units = "in", res = 300)
pheatmap(avg_exp_sub, scale = "row")
dev.off()

# =========================
# DotPlot & Signatures (PNG)
# =========================
genes_alvos <- c("CD3D", "CD3E", "MS4A1", "NKG7", "LYZ", "EPCAM", "KRT19", "COL1A1", "COL1A2", "PECAM1")
genes_present <- intersect(genes_alvos, rownames(obj))
if (length(genes_present) > 0) {
  png(file.path(outdir, "DotPlot_Genes_Alvos.png"), width = 10, height = 6, units = "in", res = 300)
  print(DotPlot(obj, features = genes_present) + RotatedAxis())
  dev.off()
}

assinatura <- list(c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4"))
assinatura_presente <- intersect(assinatura[[1]], rownames(obj))
if (length(assinatura_presente) >= 2) {
  obj <- AddModuleScore(obj, features = list(assinatura_presente), name = "Score_Assinatura")
  png(file.path(outdir, "Score_Assinatura.png"), width = 8, height = 6, units = "in", res = 300)
  print(FeaturePlot(obj, features = "Score_Assinatura1", reduction = "umap"))
  dev.off()
}

# =========================
# SingleR annotation with HPCA (Análise por Cluster - Fix)
# =========================
message("Running SingleR annotation...")
hpca.se <- celldex::HumanPrimaryCellAtlasData()

obj <- JoinLayers(obj)
obj_for_singler <- obj[["RNA"]]$data

pred.hpca <- SingleR(
  test = obj_for_singler,
  ref = hpca.se,
  labels = hpca.se$label.main,
  clusters = obj$seurat_clusters 
)

container_labels <- pred.hpca$labels[match(obj$seurat_clusters, rownames(pred.hpca))]
obj$HPCA_label <- container_labels

png(file.path(outdir, "UMAP_HPCA.png"), width = 9, height = 7, units = "in", res = 300)
print(DimPlot(obj, reduction = "umap", group.by = "HPCA_label", label = TRUE, repel = TRUE))
dev.off()

gc()
saveRDS(obj, file.path(outdir, "GSE236581_Final.rds"))
message("Analysis completed successfully.")
