library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

options(stringsAsFactors = FALSE)
set.seed(42)

data_dir <- Sys.getenv("KIDNEY_SCRNA_DATA_DIR", unset = "data")
results_dir <- Sys.getenv("KIDNEY_SCRNA_RESULTS_DIR", unset = "results")
cluster_resolution <- as.numeric(
  Sys.getenv("KIDNEY_SCRNA_CLUSTER_RESOLUTION", unset = "0.25")
)
if (is.na(cluster_resolution) || cluster_resolution <= 0) {
  stop("KIDNEY_SCRNA_CLUSTER_RESOLUTION must be a positive number.")
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- c("GSM4145204", "GSM4145205", "GSM4145206")
minimum_cells <- c(8, 6, 10)

seurat_objects <- Map(
  function(sample_id, min_cells) {
    counts <- Seurat::Read10X(data.dir = file.path(data_dir, sample_id))
    Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = min_cells,
      min.features = 200
    )
  },
  sample_ids,
  minimum_cells
)
names(seurat_objects) <- sample_ids

kidney <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = sample_ids,
  project = "human_kidney"
)
rm(seurat_objects)

kidney[["percent.mt"]] <- Seurat::PercentageFeatureSet(
  kidney,
  pattern = "^MT-"
)

qc_violin <- Seurat::VlnPlot(
  kidney,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggplot2::ggsave(
  file.path(results_dir, "qc_violin.pdf"),
  qc_violin,
  width = 12,
  height = 4
)

qc_scatter <- Seurat::FeatureScatter(
  kidney,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
) +
  Seurat::FeatureScatter(
    kidney,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
)
ggplot2::ggsave(
  file.path(results_dir, "qc_scatter.pdf"),
  qc_scatter,
  width = 10,
  height = 5
)

kidney <- subset(
  kidney,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 30
)
kidney <- Seurat::NormalizeData(kidney)
kidney <- Seurat::FindVariableFeatures(
  kidney,
  selection.method = "vst",
  nfeatures = 2000
)

variable_feature_plot <- Seurat::VariableFeaturePlot(kidney)
top_features <- head(Seurat::VariableFeatures(kidney), 10)
variable_feature_plot <- Seurat::LabelPoints(
  plot = variable_feature_plot,
  points = top_features,
  repel = TRUE
)
ggplot2::ggsave(
  file.path(results_dir, "variable_features.pdf"),
  variable_feature_plot,
  width = 8,
  height = 6
)

cell_cycle_genes <- Seurat::cc.genes.updated.2019
kidney <- Seurat::CellCycleScoring(
  kidney,
  s.features = cell_cycle_genes$s.genes,
  g2m.features = cell_cycle_genes$g2m.genes,
  set.ident = FALSE
)
kidney <- Seurat::ScaleData(
  kidney,
  vars.to.regress = c("S.Score", "G2M.Score")
)
kidney <- Seurat::RunPCA(
  kidney,
  features = Seurat::VariableFeatures(kidney),
  npcs = 20,
  verbose = FALSE
)

kidney <- harmony::RunHarmony(
  kidney,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:20,
  plot_convergence = FALSE,
  verbose = TRUE
)
kidney <- Seurat::RunUMAP(
  kidney,
  reduction = "harmony",
  dims = 1:20
)
kidney <- Seurat::FindNeighbors(
  kidney,
  reduction = "harmony",
  dims = 1:20
)
kidney <- Seurat::FindClusters(
  kidney,
  resolution = cluster_resolution
)

umap_samples <- Seurat::DimPlot(
  kidney,
  reduction = "umap",
  group.by = "orig.ident",
  split.by = "orig.ident",
  pt.size = 0.1
)
umap_phase <- Seurat::DimPlot(
  kidney,
  reduction = "umap",
  group.by = "Phase",
  pt.size = 0.1
)
umap_clusters <- Seurat::DimPlot(
  kidney,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.1
)
ggplot2::ggsave(
  file.path(results_dir, "umap_by_sample.pdf"),
  umap_samples,
  width = 12,
  height = 5
)
ggplot2::ggsave(
  file.path(results_dir, "umap_phase.pdf"),
  umap_phase,
  width = 7,
  height = 5
)
ggplot2::ggsave(
  file.path(results_dir, "umap_clusters.pdf"),
  umap_clusters,
  width = 7,
  height = 5
)

markers <- Seurat::FindAllMarkers(
  kidney,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.table(
  markers,
  file = file.path(results_dir, "cluster_markers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

marker_panel <- c(
  "SLC13A3", "SLC34A1", "GPX3", "DCXR", "SLC17A3", "SLC22A8",
  "SLC22A7", "GNLY", "NKG7", "CD3D", "CD3E", "LYZ", "CD14",
  "KRT8", "KRT18", "CD24", "VCAM1", "UMOD", "DEFB1", "CLDN8",
  "AQP2", "CD79A", "CD79B", "ATP6V1G3", "ATP6V0D2", "TMEM213"
)
marker_panel <- intersect(marker_panel, rownames(kidney))
if (length(marker_panel) > 0) {
  marker_heatmap <- Seurat::DoHeatmap(kidney, features = marker_panel)
  ggplot2::ggsave(
    file.path(results_dir, "marker_heatmap.pdf"),
    marker_heatmap,
    width = 12,
    height = 10
  )
}

kidney <- Seurat::RunTSNE(
  kidney,
  reduction = "harmony",
  dims = 1:20
)
tsne_samples <- Seurat::DimPlot(
  kidney,
  reduction = "tsne",
  group.by = "orig.ident",
  split.by = "orig.ident"
)
tsne_phase <- Seurat::DimPlot(
  kidney,
  reduction = "tsne",
  group.by = "Phase"
)
ggplot2::ggsave(
  file.path(results_dir, "tsne_by_sample.pdf"),
  tsne_samples,
  width = 12,
  height = 5
)
ggplot2::ggsave(
  file.path(results_dir, "tsne_phase.pdf"),
  tsne_phase,
  width = 7,
  height = 5
)

saveRDS(kidney, file.path(results_dir, "integrated_kidney_seurat.rds"))

proximal_tubule_clusters <- intersect(
  c("0", "1", "2"),
  levels(Seurat::Idents(kidney))
)
if (length(proximal_tubule_clusters) > 0) {
  proximal_tubule <- subset(
    kidney,
    idents = proximal_tubule_clusters
  )
  saveRDS(
    proximal_tubule,
    file.path(results_dir, "proximal_tubule_subset.rds")
  )
}
