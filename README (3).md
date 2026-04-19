# 🧬 Single-Cell RNA-seq Tutorial Collection

> A structured reference guide covering three foundational tutorials for scRNA-seq preprocessing, data structures, and downstream analysis using Galaxy, AnnData, and Scanpy/scverse.

---

## 📚 Table of Contents

- [Overview](#overview)
- [Tutorial 1 — 10X scRNA Preprocessing (Galaxy)](#tutorial-1--pre-processing-of-10x-single-cell-rna-datasets-galaxy)
- [Tutorial 2 — Getting Started with AnnData](#tutorial-2--getting-started-with-anndata)
- [Tutorial 3 — Basic scRNA Analysis (scverse / Scanpy)](#tutorial-3--basic-scrna-tutorial-scverse--scanpy)
- [Tool Summary Across Tutorials](#tool-summary-across-tutorials)
- [Recommended Learning Order](#recommended-learning-order)

---

## Overview

| # | Tutorial | Platform | Dataset | Difficulty |
|---|----------|----------|---------|------------|
| 1 | [Pre-processing of 10X scRNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html) | Galaxy (web-based, no coding) | 1k PBMCs from 10x Genomics (v3) | Beginner |
| 2 | [Getting Started with AnnData](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) | Python / Jupyter Notebook | Simulated (100 cells × 2000 genes) | Beginner–Intermediate |
| 3 | [Basic scRNA Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html) | Python / Jupyter (Scanpy + scverse) | Multi-sample PBMC dataset | Intermediate |

---
---

# Tutorial 1 — Pre-processing of 10X Single-Cell RNA Datasets (Galaxy)

🔗 **Link:** https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html

**Authors:** Wendi Bacon, Pavankumar Videm, Mehmet Tekman, Hans-Rudolf Hotz, Daniel Blankenberg
**Last Updated:** April 29, 2025 | **License:** CC-BY-4.0

### 🎯 Objective
Learn to preprocess 10x Genomics scRNA-seq FASTQ files into a filtered count matrix using the Galaxy platform — no programming required.

### 🧾 Dataset
- **Source:** 1k PBMCs from a Healthy Donor (10x Genomics, v3 chemistry)
- **Tutorial subset:** ~300 cells (sub-sampled for training speed)
- **Input files (6 FASTQ):**

| File | Lane | Read Type | Content |
|------|------|-----------|---------|
| `L001_R1_001.fastq.gz` | L001 | R1 | Cell Barcodes |
| `L001_R2_001.fastq.gz` | L001 | R2 | cDNA Sequences |
| `L001_I1_001.fastq.gz` | L001 | I1 | Illumina Lane Info |
| `L002_R1_001.fastq.gz` | L002 | R1 | Cell Barcodes |
| `L002_R2_001.fastq.gz` | L002 | R2 | cDNA Sequences |
| `L002_I1_001.fastq.gz` | L002 | I1 | Illumina Lane Info |

> ⚠️ RNA STARsolo does **not** need the I1 files. They are used by Cell Ranger only.

---

### 🔬 Step-by-Step Workflow

#### Step 1 — Data Upload to Galaxy

| Field | Value |
|-------|-------|
| **Action** | Upload 6 FASTQ files + GTF annotation file + barcode whitelist |
| **Input** | FASTQ (R1, R2) + `Homo_sapiens.GRCh37.75.gtf` + `3M-february-2018.txt.gz` |
| **Output** | Files available in Galaxy History panel |

---

#### Step 2 — Mapping & Quantification with RNA STARsolo

| Field | Value |
|-------|-------|
| **Tool** | `RNA STARsolo` (Galaxy version 2.7.10b+galaxy3) |
| **Input** | R1 FASTQ (barcodes): L001 + L002 · R2 FASTQ (cDNA): L001 + L002 · Reference genome · GTF · Barcode whitelist |
| **Key Parameters** | Reference: hg19 (chrX for training) · Junction overhang: 100bp · SC type: Drop-seq / 10X Chromium |
| **Output** | Count matrix · Feature Statistic Summaries · BAM file · Log files |
| **Result** | ~5200 raw barcodes detected; 272 pass default quality threshold |

> 💡 STARsolo is a drop-in replacement for Cell Ranger — much faster and no complex configuration needed.

---

#### Step 3 — Quality Control Inspection

| Field | Value |
|-------|-------|
| **Tool** | Galaxy eye viewer on STARsolo Feature Statistic Summaries |
| **Input** | STARsolo log + Feature Statistic Summaries file |
| **Output** | QC report showing `yesCellBarcodes` count |
| **Result** | Confirms 272 high-quality barcodes above threshold (from 5200 total) |

---

#### Step 4 — Barcode Ranking (Knee Plot)

| Field | Value |
|-------|-------|
| **Tool** | `DropletUtils` |
| **Input** | STARsolo count matrix |
| **Output** | Knee/inflection plot; ranked barcode table |
| **Result** | Visual "knee" separating real cells from empty droplets |

---

#### Step 5 — Custom Filtering with DropletUtils

| Field | Value |
|-------|-------|
| **Tool** | `DropletUtils` (custom threshold mode) |
| **Input** | Count matrix + thresholds from knee plot |
| **Output** | Filtered count matrix |
| **Result** | 282 high-quality cells (vs. 272 from default method) |

---

### 📊 Workflow Summary

```
FASTQ Files (R1 + R2, 2 lanes)
        │
        ▼
  [RNA STARsolo]  ← Reference Genome + GTF + Barcode Whitelist
        │
        ├──► BAM File
        ├──► QC Log / Feature Summary
        └──► Raw Count Matrix
                  │
                  ▼
          [DropletUtils]
                  │
                  ├──► Knee Plot (barcode ranking)
                  └──► Filtered Count Matrix (282 cells)
```

---

### 📦 Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| RNA STARsolo | 2.7.10b+galaxy3 | Alignment + UMI quantification |
| DropletUtils | Galaxy latest | Barcode filtering + knee plot |
| Galaxy Log Viewer | — | QC visualization |

---
---

# Tutorial 2 — Getting Started with AnnData

🔗 **Link:** https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

**Authors:** Adam Gayoso, Alex Wolf
**Install:** `pip install anndata` or `conda install anndata -c conda-forge`

---

### 🎯 What Is AnnData and Why Does It Matter?

AnnData (Annotated Data) is the **standard data container** for single-cell analysis in Python. Almost every major Python-based scRNA-seq tool — Scanpy, scVI, CellTypist, Decoupler — stores and reads data as AnnData objects.

The problem it solves: in scRNA-seq, you don't just have a matrix of numbers. For each cell you might have metadata like cell type, donor, batch, or QC scores. For each gene you might have chromosome location, gene ID, or variability flags. You also need to store computed results like UMAP embeddings, PCA coordinates, and nearest-neighbor graphs — all linked consistently to the same cells and genes. AnnData keeps all of this in one organized object so nothing gets out of sync.

Think of AnnData as a spreadsheet where:
- **Rows** = cells (observations)
- **Columns** = genes (variables)
- Every annotation, embedding, graph, and result stays attached to those same rows/columns

---

### 🏗️ AnnData Object Structure

```
AnnData  (n_obs × n_vars = 100 × 2000)
│
├── .X           ← The core count matrix (cells × genes), stored as sparse matrix
├── .obs         ← Cell-level metadata table (Pandas DataFrame)
├── .var         ← Gene-level metadata table (Pandas DataFrame)
├── .obsm        ← Multi-dimensional cell data (e.g., UMAP coords, PCA scores)
├── .varm        ← Multi-dimensional gene data (e.g., PCA loadings)
├── .obsp        ← Cell-cell pairwise data (e.g., KNN distance matrix)
├── .varp        ← Gene-gene pairwise data (e.g., co-expression graph)
├── .layers      ← Alternative versions of .X (e.g., log-normalized, raw)
└── .uns         ← Unstructured metadata (dicts, lists, color maps, etc.)
```

---

### 🧾 Dataset Used in This Tutorial

| Property | Value |
|----------|-------|
| Type | Simulated / synthetic (randomly generated) |
| Dimensions | 100 cells × 2000 genes |
| Distribution | Poisson(λ=1) — mimics sparse count data typical in scRNA-seq |
| Why synthetic? | Lets you focus on the data structure itself, not the biology |

---

### 🔬 Step-by-Step Walkthrough

---

#### Step 1 — Installation & Imports

```python
pip install anndata numpy pandas scipy
```

```python
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
```

**Why these libraries?**
`numpy` and `pandas` handle numerical and tabular operations. `scipy.sparse.csr_matrix` creates a sparse matrix — this is important because scRNA-seq data is extremely sparse (most genes have zero expression in most cells), so storing it as a dense matrix would waste enormous memory. A sparse format only stores the non-zero values.

---

#### Step 2 — Creating the Core Count Matrix and AnnData Object

```python
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
# AnnData object with n_obs × n_vars = 100 × 2000
```

**What's happening here:**
`np.random.poisson(1, size=(100, 2000))` generates a 100×2000 matrix of random integer counts drawn from a Poisson distribution — this mimics real gene expression counts, which are naturally count data and highly sparse. `csr_matrix(...)` converts it to Compressed Sparse Row format for memory efficiency. `ad.AnnData(counts)` wraps this matrix into an AnnData object, accessible as `adata.X`.

At this point the object has no labels — cells are just numbered 0–99 and genes 0–1999. You'd always want to name them next.

---

#### Step 3 — Naming Cells and Genes (`obs_names` / `var_names`)

```python
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
```

**Why this matters:**
Every cell and gene must have a **unique string identifier**. These names serve as the shared index across all annotations, embeddings, and graphs inside the object. When you compute a UMAP, filter cells, or look up gene metadata later, AnnData uses these names to keep everything aligned. Functions in the scverse ecosystem raise warnings if `obs_names` or `var_names` are not unique.

In real data, `obs_names` are typically **cell barcodes** (e.g., `AAACATACAACCAC-1`) and `var_names` are **gene IDs or gene symbols** (e.g., `BRCA1`, `ENSG00000012048`).

---

#### Step 4 — Adding Cell-Level Metadata (`.obs`)

```python
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)
```

**What's happening and why:**
`.obs` is a Pandas DataFrame where each row is a cell. You can add any number of columns to it — cell type, donor ID, batch, sequencing depth, QC pass/fail flags, etc. Adding metadata here keeps it permanently linked to the cells it describes, so if you later filter or subset cells, the metadata follows automatically.

Using `pd.Categorical` is recommended for columns with repeated string values (like cell types) because it saves memory by storing the unique categories once and referring to them by integer codes internally.

In a real scRNA workflow, `.obs` grows throughout the analysis — first you add QC metrics, then doublet scores, then cluster assignments, then final cell type annotations.

```
         cell_type
Cell_0    Monocyte
Cell_1         T
Cell_2         B
```

---

#### Step 5 — Adding Gene-Level Metadata (`.var`)

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # real-world example
```

**What's happening and why:**
`.var` works exactly like `.obs` but for genes. In real data you store things like: whether a gene is highly variable, whether it's a mitochondrial gene, its chromosome location, or alternative gene symbols. For example, Scanpy's `pp.highly_variable_genes()` adds a column `highly_variable = True/False` to `.var`, and PCA automatically uses this column to restrict computation to only those genes — so the metadata in `.var` directly drives downstream analysis steps.

---

#### Step 6 — Storing Multi-Dimensional Embeddings (`.obsm` and `.varm`)

```python
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))
```

**Why `.obsm` instead of `.obs`?**
`.obs` can only hold one value per cell per column (a scalar or a string). But dimensionality reduction produces a whole vector per cell — UMAP gives 2 coordinates, PCA might give 50 components. You cannot put a 50-dimensional vector into a single DataFrame cell.

`.obsm` (observation matrix) is a dictionary of arrays where each value has shape `(n_cells, k)`. The key name is arbitrary — by convention:
- `"X_umap"` → 2D UMAP coordinates
- `"X_pca"` → PCA embedding (typically 50 dimensions)
- `"X_tsne"` → t-SNE coordinates

`.varm` works the same way for genes — commonly used to store PCA loadings, which describe how much each gene contributes to each principal component.

---

#### Step 7 — Storing Pairwise Relationships (`.obsp` and `.varp`)

```python
# .obsp holds matrices of shape (n_cells × n_cells)
adata.obsp["distances"] = sparse_distance_matrix  # shape: (100, 100)
```

**Why is this needed?**
Many scRNA algorithms operate on graphs rather than on gene expression directly. When Scanpy's `pp.neighbors()` runs, it computes a distance matrix and connectivity matrix between all pairs of cells and stores them in `.obsp["distances"]` and `.obsp["connectivities"]`. Leiden clustering and UMAP then use this graph — so `.obsp` is a critical intermediate storage location in any standard pipeline.

`.varp` does the equivalent for gene-gene relationships, used for gene co-expression or regulatory network analysis.

---

#### Step 8 — Storing Alternative Data in Layers (`.layers`)

```python
adata.layers["log_transformed"] = np.log1p(adata.X)
```

**What is a layer and why use it?**
Layers let you store multiple versions of your count matrix at the same time. During a typical pipeline, you may want to preserve the raw integer counts while also having normalized and log-transformed versions. Instead of creating separate objects (which risk getting out of sync), AnnData layers store all versions together, each with the same shape and aligned to the same cells and genes.

In practice:
- `adata.layers["counts"]` → original raw integer counts (preserved)
- `adata.X` → log-normalized counts (actively used by analysis functions)
- `adata.layers["scaled"]` → zero-mean, unit-variance matrix (used by PCA)

You can retrieve any layer as a Pandas DataFrame: `adata.to_df(layer="log_transformed")`.

---

#### Step 9 — Storing Unstructured Metadata (`.uns`)

```python
adata.uns["random"] = {"arbitrary_key": "arbitrary_value"}
```

**What goes in `.uns`?**
`.uns` (unstructured) is a free-form dictionary for anything that doesn't fit neatly into the other slots. Typical uses include: analysis parameters (e.g., number of neighbors used for KNN graph), color palettes for cell type plotting, differential expression results from `rank_genes_groups`, and general experiment metadata like species, tissue, or protocol information.

---

#### Step 10 — Subsetting the Object

```python
# Subset by cell name
subset = adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]

# Subset by boolean condition
t_cells = adata[adata.obs["cell_type"] == "T"]
```

**Views vs. Copies — an important distinction:**
Subsetting returns a **view**, not a copy. A view is a window into the original object — no data is duplicated in memory. This is efficient when you just want to inspect a subset. However, if you want an independent object that you can modify without affecting the original, you must explicitly call `.copy()`:
```python
t_cells = adata[adata.obs["cell_type"] == "T"].copy()
```
This distinction matters a lot when running filters in a pipeline — always `.copy()` after a major filter step.

---

#### Step 11 — Saving and Loading (`.h5ad` format)

```python
adata.write_h5ad("my_data.h5ad")     # Save entire object to disk
adata = ad.read_h5ad("my_data.h5ad") # Load it back
```

**What is `.h5ad`?**
AnnData's native file format is `.h5ad`, an HDF5-based binary file. HDF5 supports hierarchical storage of arrays, DataFrames, and dictionaries — a perfect fit for the AnnData structure. Saving to `.h5ad` preserves ALL slots (`.X`, `.obs`, `.var`, `.obsm`, `.obsp`, `.layers`, `.uns`) in one portable file, making it easy to share data and resume analysis without rerunning computations.

---

### 📊 AnnData Attribute Reference Table

| Attribute | Shape | Type | What It Stores | Real Example |
|-----------|-------|------|---------------|--------------|
| `.X` | (n_obs × n_vars) | Sparse/dense matrix | Core expression counts | Raw UMI counts |
| `.obs` | (n_obs × cols) | Pandas DataFrame | Cell metadata | `cell_type`, `batch`, `n_genes`, `pct_mt` |
| `.var` | (n_vars × cols) | Pandas DataFrame | Gene metadata | `highly_variable`, `mean`, `chromosome` |
| `.obsm` | (n_obs × k) per key | Dict of arrays | Multi-dim cell data | `X_umap`, `X_pca`, `X_tsne` |
| `.varm` | (n_vars × k) per key | Dict of arrays | Multi-dim gene data | `PCs` (PCA loadings) |
| `.obsp` | (n_obs × n_obs) per key | Dict of sparse matrices | Cell-cell relationships | `distances`, `connectivities` |
| `.varp` | (n_vars × n_vars) per key | Dict of sparse matrices | Gene-gene relationships | Gene regulatory network |
| `.layers` | (n_obs × n_vars) per key | Dict of matrices | Alternative count matrices | `counts`, `log1p`, `scaled` |
| `.uns` | Unstructured | Dict | Free-form metadata | Analysis params, color maps, DE results |

---
---

# Tutorial 3 — Basic scRNA Tutorial (scverse / Scanpy)

🔗 **Link:** https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html

**Framework:** Scanpy + scverse ecosystem
**Key Packages:** `scanpy`, `anndata`, `decoupler`, `celltypist`

---

### 🎯 What Does This Tutorial Cover?

This tutorial walks through a complete, real-world scRNA-seq analysis pipeline — from a raw count matrix all the way to biologically annotated cell clusters. Unlike Tutorial 1 (which produces the count matrix) and Tutorial 2 (which teaches the data structure), this tutorial focuses on **what you actually do with the data after you have it** — cleaning it, reducing its dimensions, finding groups of similar cells, and figuring out which cell types those groups represent.

It follows current best practices in the scverse ecosystem and uses multi-sample PBMC data, so it also addresses batch handling.

---

### 🧾 Dataset

| Property | Value |
|----------|-------|
| Cell type | PBMCs (Peripheral Blood Mononuclear Cells) from healthy donors |
| Format | Pre-loaded `.h5ad` file |
| Samples | Multiple batches merged into one AnnData object |
| Pre-existing annotations | Sample ID, basic cell metadata |

---

### 🔬 Step-by-Step Walkthrough

---

#### Step 1 — Installation & Setup

```python
pip install scanpy celltypist decoupler
```

```python
import scanpy as sc
import decoupler as dc
import celltypist
```

**What these packages do:**
- `scanpy` — the main toolkit: preprocessing, dimensionality reduction, clustering, visualization — all operating on AnnData objects
- `decoupler` — computes transcription factor and pathway activity scores from gene expression data, adding biological interpretation
- `celltypist` — automated cell type annotation using models pretrained on large annotated single-cell datasets

---

#### Step 2 — Load the Data

```python
adata = sc.read_h5ad("pbmc_data.h5ad")
print(adata)
```

The `.h5ad` file already contains the gene × cell count matrix plus basic metadata (sample IDs, etc.). After loading, you inspect the object to understand: how many cells and genes it contains, what columns already exist in `.obs`, and whether the data looks reasonable before any processing.

---

#### Step 3 — Quality Control (QC)

```python
adata.var["mt"]   = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"]   = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)
```

**Why do QC?**
Real scRNA-seq data always contains low-quality cells that must be removed before any analysis — including them would create false clusters or corrupt genuine ones. There are three main sources of bad data:

1. **Dead or dying cells** — damaged cells lose cytoplasmic RNA but retain mitochondrial RNA, so a high fraction of mitochondrial gene reads (`pct_counts_mt`) signals cell death. These cells are typically filtered out above a threshold of ~20%.

2. **Empty droplets** — sometimes a droplet captures no cell at all but still gets a barcode. These produce cells with very few detected genes. Filtering cells with fewer than ~200 genes removes them.

3. **Doublets** — two cells captured in one droplet, appearing as a single artificially large cell with unusually high gene counts. Handled more precisely in Step 5.

The first part of QC is not filtering yet — it is **computing and labeling** the metrics so you can visualize them and decide appropriate thresholds.

| QC Metric | What It Measures | Action |
|-----------|-----------------|--------|
| `n_genes_by_counts` | Genes detected per cell | Too low = empty droplet · Too high = doublet |
| `total_counts` | Total UMI counts per cell | Correlates with depth; very high suggests doublet |
| `pct_counts_mt` | % reads from mitochondrial genes | High (>20%) = dying/damaged cell → remove |
| `pct_counts_ribo` | % reads from ribosomal genes | Context-dependent |
| `pct_counts_hb` | % reads from hemoglobin genes | Should be near zero in PBMCs |

---

#### Step 4 — Filtering Cells and Genes

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
```

**What each filter does and why:**

- `min_genes=200`: Removes cells that express fewer than 200 genes. No real cell should express fewer than a few hundred genes — these are empty droplets or failed captures. Keeping them would pollute the analysis with noise.

- `min_cells=3`: Removes genes that are detected in fewer than 3 cells across the entire dataset. These genes are far too rare to contribute to clustering or to be statistically meaningful. Removing them reduces the matrix size and speeds up computation without any loss of information.

- `pct_counts_mt < 20`: Removes cells where more than 20% of reads come from mitochondrial genes. These are likely damaged or dying cells.

The exact thresholds vary by tissue, species, and protocol — you should always visualize the QC metrics first and choose cutoffs that make biological sense for your specific dataset.

---

#### Step 5 — Doublet Detection with Scrublet

```python
sc.pp.scrublet(adata)
```

**What is a doublet?**
A doublet occurs when two cells are co-captured in one droplet and sequenced as a single unit. They appear as cells with roughly double the gene count of typical cells and can create spurious clusters or distort real ones in downstream analysis — so detecting and removing them is important.

**How Scrublet works:**
Scrublet simulates artificial doublets by randomly combining pairs of real cells from your data (since a real doublet is essentially two cells merged together). It then trains a nearest-neighbor classifier to score how closely each observed cell resembles these simulated doublets.

The output — a `doublet_score` and binary `predicted_doublet` flag — is added to `adata.obs`. You can then either:
- Remove predicted doublets immediately: `adata = adata[~adata.obs["predicted_doublet"]].copy()`
- Or use the scores later after clustering to identify and remove entire doublet-enriched clusters

---

#### Step 6 — Normalization

```python
sc.pp.normalize_total(adata, target_sum=1e4)  # per-cell normalization
sc.pp.log1p(adata)                             # log(x + 1) transformation
adata.layers["log1p_norm"] = adata.X.copy()   # save this version as a layer
```

**Why normalize?**
Raw count data has a major technical problem: different cells were sequenced to different depths. One cell might have 5,000 total UMI counts and another might have 20,000 — not because the second cell actually expresses more, but because it happened to be sequenced more deeply. If you don't correct for this, highly sequenced cells will falsely appear more similar to each other than to shallowly sequenced cells of the same type.

**Step 1 — Count depth normalization (`normalize_total`):** Rescales each cell's counts so that every cell sums to the same total (10,000 counts). After this, a gene with count 50 means the same proportion of expression in every cell, making cells comparable.

**Step 2 — Log1p transformation:** After normalization, gene expression is still heavily right-skewed — a few highly expressed genes have values hundreds of times larger than most genes. Log-transforming with `log(x + 1)` (the +1 prevents log(0)) compresses this scale, making expression values more normally distributed and better suited for PCA and other linear methods that assume approximate normality.

After saving the normalized version as `layers["log1p_norm"]`, the raw integer counts should also be preserved (typically in `layers["counts"]`) for methods that require them.

---

#### Step 7 — Feature Selection (Highly Variable Genes)

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
```

**Why select features?**
Human cells express roughly 20,000 genes, but the vast majority either are not expressed at all or are expressed at a constant level across all cells (housekeeping genes like ribosomal proteins). These uninformative genes add noise without contributing any information about cell type differences.

Highly Variable Genes (HVGs) are genes whose expression level differs significantly across cells — they are the genes that actually distinguish one cell type from another. By selecting only the top 2,000 most variable genes, you reduce the working dimensionality from ~20,000 to 2,000, remove noise, and make downstream computations (PCA, UMAP, clustering) much faster and more meaningful.

The `batch_key="sample"` argument is important for multi-sample data: it selects genes that are variable *within* each sample before combining results across samples. This prevents the selection from being dominated by genes that merely differ between batches rather than between real cell types.

After this step, `adata.var["highly_variable"]` is a boolean column. All subsequent steps (scaling, PCA) automatically operate only on these 2,000 genes.

---

#### Step 8 — Scaling

```python
sc.pp.scale(adata, max_value=10)
```

**Why scale?**
Even after normalization and log transformation, different genes still have very different mean expression levels and variances. PCA is sensitive to absolute scale — a gene with large variance would dominate the principal components simply because its values are larger, not because it is biologically more important.

Scaling transforms each gene to have zero mean and unit variance across all cells. After scaling, each gene contributes equally to PCA regardless of its original expression level. The `max_value=10` clips extreme outliers to prevent them from distorting the scaled space.

---

#### Step 9 — Dimensionality Reduction with PCA

```python
sc.pp.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50)
```

**What is PCA doing here?**
Even with 2,000 genes, the data is still very high-dimensional, and much of that variation is redundant (correlated genes carry similar information). Principal Component Analysis (PCA) finds the directions in gene-expression space that explain the most variation and represents each cell as a point in this lower-dimensional space.

The result is stored in `adata.obsm["X_pca"]` — each cell now has a vector of ~50 PCA scores instead of 2,000 gene values. This compact representation captures most of the biological variation while dramatically reducing noise and computation time.

The **variance ratio plot** (elbow plot) shows how much variance each principal component explains. You look for the "elbow" where the curve flattens — this tells you how many PCs carry real biological signal (typically 15–40 PCs). Using too few PCs loses signal; using too many adds noise back in.

---

#### Step 10 — Building the Neighborhood Graph (KNN Graph)

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
```

**What is the neighborhood graph?**
This is the most critical structural step before clustering and UMAP. The concept is simple: cells that are similar in PCA space should be connected. For each cell, Scanpy finds its `n_neighbors` (15 here) nearest neighbors in PCA space and connects them with edges, building a graph where nodes are cells and edges connect similar cells.

This graph is the foundation for all downstream steps — both UMAP and Leiden clustering operate on this graph. The results are stored in `adata.obsp["distances"]` and `adata.obsp["connectivities"]`.

The `n_pcs` parameter tells Scanpy how many principal components to use when computing distances. This should match what you determined from the elbow plot in Step 9.

---

#### Step 11 — UMAP Visualization

```python
sc.tl.umap(adata)
sc.pl.umap(adata, color="sample")
```

**What is UMAP?**
UMAP (Uniform Manifold Approximation and Projection) takes the neighborhood graph and projects cells into 2D in a way that preserves local structure — cells that are neighbors in high-dimensional PCA space stay close together in 2D. The result is a scatter plot where groups of similar cells appear as visual islands or clusters.

UMAP is used for **visualization only** — the actual clustering (Step 12) is done on the graph, not the UMAP coordinates. The 2D layout is not quantitatively reliable for measuring distances between cells (two islands might be far apart on UMAP even if they share many marker genes). However, UMAP is invaluable for checking that your preprocessing worked and for visualizing cluster assignments and gene expression patterns.

The 2D coordinates are stored in `adata.obsm["X_umap"]`.

---

#### Step 12 — Clustering (Leiden Algorithm)

```python
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res0_5")
sc.pl.umap(adata, color="leiden_res0_5", legend_loc="on data")
```

**What is Leiden clustering?**
Leiden is a graph community detection algorithm that partitions the neighborhood graph into groups of highly connected cells. Unlike k-means, you do not need to specify the number of clusters in advance. Instead, you control granularity with the `resolution` parameter:

- **Low resolution** (e.g., 0.1–0.3) → fewer, larger clusters capturing major cell lineages
- **High resolution** (e.g., 0.8–1.5) → more, smaller clusters resolving fine cell subtypes

The cluster labels are added as a column in `adata.obs`. It is common practice to run Leiden at multiple resolutions and compare before settling on one that matches your biological question.

---

#### Step 13 — Identifying Marker Genes

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_res0_5", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
```

**What are marker genes?**
Now that you have clusters, you need to understand what biological cell type each cluster represents. Marker genes are genes that are significantly more highly expressed in one cluster compared to all others — they are the molecular "signature" of each cluster.

`rank_genes_groups` performs a Wilcoxon rank-sum test (recommended for scRNA-seq data) for each gene in each cluster versus all other clusters. Results are stored in `adata.uns["rank_genes_groups"]` as ranked lists of differentially expressed genes per cluster.

The **dotplot** visualization provides a compact summary:
- **Dot size** = fraction of cells in the cluster that express the gene (detection rate)
- **Dot color** = mean expression level in expressing cells

Together these help you identify the distinguishing genes per cluster and look them up in known cell type databases.

---

#### Step 14 — Cell Type Annotation

Once marker genes are identified, you assign biological identities to the clusters. Three methods are shown in this tutorial:

**Method A — Manual annotation using known marker genes:**
```python
marker_genes = {
    "B cells": ["MS4A1", "CD79A"],
    "T cells": ["CD3D", "CD3E"],
    "Monocytes": ["LYZ", "CD14"],
    "NK cells": ["NKG7", "GNLY"]
}
sc.pl.umap(adata, color=["MS4A1", "CD3D", "LYZ", "NKG7"])
```
You overlay the expression of known marker genes on the UMAP. Clusters where `MS4A1` is highly expressed are B cells, clusters where `CD3D` is high are T cells, and so on. This approach requires prior biological knowledge of what cell types to expect in the tissue.

**Method B — Automated annotation with CellTypist:**
```python
predictions = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
adata.obs["celltypist_cell_type"] = predictions.predicted_labels
```
CellTypist uses a logistic regression model trained on thousands of annotated single-cell immune datasets to predict cell type labels for each cell automatically. The `majority_voting=True` option takes the most common prediction within each Leiden cluster as the consensus label, smoothing out individual cell noise. This is useful as a first-pass annotation especially when you are less familiar with the expected cell types.

**Method C — Transcription factor activity scoring with Decoupler:**
```python
dc.run_mlm(mat=adata, net=net, source="source", target="target", weight="weight")
acts = dc.pp.get_obsm(adata, "score_mlm")
sc.pl.umap(acts, color=["majority_voting", "B cells", "T cells", "Monocytes"])
```
Instead of looking at raw gene expression, Decoupler infers the activity of transcription factors (TFs) or biological pathways from the expression patterns of their target genes. This adds a mechanistic layer of interpretation — rather than just observing that cluster 3 expresses MS4A1, you can confirm that it has high activity of the B cell master regulator PAX5, strengthening the annotation. Decoupler results also help discover regulatory biology beyond just assigning cell type names.

The three methods are complementary: use CellTypist for a quick automated first pass, confirm with manual markers on the UMAP, and use Decoupler to gain mechanistic insight.

---

### 📊 Full Pipeline Summary Table

| Step | Function | What It Does | Result Stored In |
|------|----------|-------------|-----------------|
| 1. Load | `sc.read_h5ad()` | Read count matrix + metadata from disk | `adata` object |
| 2. Tag genes | `adata.var["mt"] = ...` | Label mitochondrial / ribosomal / hemoglobin genes | `adata.var` new columns |
| 3. QC metrics | `sc.pp.calculate_qc_metrics()` | Compute per-cell quality statistics | `adata.obs` new columns |
| 4. Filter | `sc.pp.filter_cells/genes()` | Remove low-quality cells and rare genes | Smaller filtered `adata` |
| 5. Doublets | `sc.pp.scrublet()` | Score each cell as likely singlet or doublet | `adata.obs["doublet_score"]` |
| 6. Normalize | `normalize_total()` + `log1p()` | Remove sequencing depth bias; stabilize variance | Updated `adata.X` |
| 7. HVG | `sc.pp.highly_variable_genes()` | Select 2000 most variable genes for analysis | `adata.var["highly_variable"]` |
| 8. Scale | `sc.pp.scale()` | Zero mean, unit variance per gene | Updated `adata.X` |
| 9. PCA | `sc.pp.pca()` | Compress 2000 genes → 50 PCA dimensions | `adata.obsm["X_pca"]` |
| 10. Neighbors | `sc.pp.neighbors()` | Build KNN graph connecting similar cells | `adata.obsp["distances/connectivities"]` |
| 11. UMAP | `sc.tl.umap()` | Project cells to 2D for visualization | `adata.obsm["X_umap"]` |
| 12. Cluster | `sc.tl.leiden()` | Partition graph into cell communities | `adata.obs["leiden_res0_5"]` |
| 13. Markers | `sc.tl.rank_genes_groups()` | Find genes distinguishing each cluster | `adata.uns["rank_genes_groups"]` |
| 14. Annotate | CellTypist / Decoupler / manual | Assign biological cell type names to clusters | `adata.obs["cell_type"]` |

---

### 📦 Packages Used

| Package | Role |
|---------|------|
| `scanpy` | Core single-cell toolkit: preprocessing, clustering, visualization |
| `anndata` | Data container (AnnData object used throughout) |
| `celltypist` | Automated cell type annotation with pretrained models |
| `decoupler` | Transcription factor and pathway activity scoring |
| `numpy` / `pandas` | Numerical computation and tabular data handling |

---
---

## Tool Summary Across Tutorials

| Tool / Package | Tutorial 1 (Galaxy) | Tutorial 2 (AnnData) | Tutorial 3 (Scanpy) |
|----------------|--------------------|-----------------------|---------------------|
| RNA STARsolo | ✅ Core aligner | ❌ | ❌ |
| DropletUtils | ✅ Barcode filtering | ❌ | ❌ |
| AnnData | Output format | ✅ Core focus | ✅ Data container |
| Scanpy | ❌ | ❌ | ✅ Core framework |
| Scrublet | ❌ | ❌ | ✅ Doublet detection |
| CellTypist | ❌ | ❌ | ✅ Cell type annotation |
| Decoupler | ❌ | ❌ | ✅ TF activity scoring |
| Galaxy Platform | ✅ | ❌ | ❌ |
| Python / Jupyter | ❌ | ✅ | ✅ |

---

## Recommended Learning Order

```
[Tutorial 1]  →  [Tutorial 2]  →  [Tutorial 3]
 Galaxy 10X        AnnData           Scanpy
 Preprocessing     Data Structure    Full Analysis
 (No coding)       (Python basics)   (Advanced Python)
```

1. **Start with Tutorial 1** if you are new to scRNA-seq — covers the biology and data flow from raw FASTQ to count matrix without any programming.
2. **Do Tutorial 2** to understand the AnnData object, which is the universal data container for all Python-based tools. You cannot use Tutorial 3 effectively without understanding `.obs`, `.obsm`, `.layers`, etc.
3. **Proceed to Tutorial 3** for a complete end-to-end Python analysis pipeline following current best practices — from loading counts through clustering and biological annotation.

---

## 📖 References

- Bacon et al. (2025). *Pre-processing of 10X Single-Cell RNA Datasets*. Galaxy Training Network. https://gxy.io/GTN:W00210
- Virshup et al. (2024). *anndata: Annotated Data*. JOSS. doi: 10.21105/joss.04371
- Wolf et al. (2018). *SCANPY: large-scale single-cell gene expression data analysis*. Genome Biology. doi: 10.1186/s13059-017-1382-0
- Tekman et al. (2020). *A single-cell RNA-sequencing training and analysis suite using the Galaxy framework*. GigaScience. doi: 10.1093/gigascience/giaa102

---

*README compiled from official tutorial documentation. For issues or updates, refer to the original tutorial links above.*
