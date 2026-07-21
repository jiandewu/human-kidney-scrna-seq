# Human Kidney scRNA-seq Analysis

A Seurat and Harmony workflow for quality control, integration, clustering,
marker discovery, and visualization of three human kidney single-cell RNA-seq
samples: GSM4145204, GSM4145205, and GSM4145206.

## Workflow

`analysis/kidney_scrna_workflow.R` performs:

1. loading and merging three 10x Genomics count matrices;
2. mitochondrial-content and feature-count quality control;
3. log normalization and variable-feature selection;
4. cell-cycle scoring and regression;
5. PCA followed by Harmony correction across sample identities;
6. UMAP, shared-nearest-neighbor clustering, and t-SNE;
7. cluster marker detection and marker-panel visualization; and
8. saving the integrated Seurat object and an exploratory proximal-tubule
   subset.

## Data policy

Count matrices, cell barcodes, sample metadata, Seurat objects, plots, and
marker tables are not tracked in this repository.

Place independently obtained 10x-format directories under `data/`:

```text
data/
|-- GSM4145204/
|   |-- barcodes.tsv.gz
|   |-- features.tsv.gz
|   `-- matrix.mtx.gz
|-- GSM4145205/
|   |-- barcodes.tsv.gz
|   |-- features.tsv.gz
|   `-- matrix.mtx.gz
`-- GSM4145206/
    |-- barcodes.tsv.gz
    |-- features.tsv.gz
    `-- matrix.mtx.gz
```

Uncompressed 10x files are also accepted by Seurat's `Read10X()`.

## Environment

```bash
conda env create -f environment.yml
conda activate human-kidney-scrna-seq
```

## Configuration

The workflow accepts the following environment variables:

```bash
export KIDNEY_SCRNA_DATA_DIR=/path/to/input
export KIDNEY_SCRNA_RESULTS_DIR=/path/to/results
export KIDNEY_SCRNA_CLUSTER_RESOLUTION=0.25
```

Defaults are `data`, `results`, and `0.25`.

## Run

From the repository root:

```bash
Rscript analysis/kidney_scrna_workflow.R
```

The workflow writes QC plots, embeddings, marker tables, and RDS objects to
the configured results directory.

## Reproducibility notes

- QC thresholds retain the original workflow values: 200–2,500 detected
  features and less than 30% mitochondrial expression.
- Harmony removes variation associated with `orig.ident`; downstream UMAP,
  neighbor, clustering, and t-SNE calculations use the Harmony embedding.
- Cluster identities remain numerical. Cell-type annotation should be
  validated from marker expression before biological interpretation.
- The optional proximal-tubule object uses clusters 0, 1, and 2 as an
  exploratory starting point. Those cluster numbers are not guaranteed to
  represent proximal-tubule cells in every run.
- Random seed 42 is used for reproducible dimensional-reduction initialization.

## Privacy

The retained history contains source code only. No count matrices, cell
barcodes, phenotype metadata, generated objects, private paths, or manuscripts
are included.
