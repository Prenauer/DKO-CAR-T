# DKO CAR-T manuscript code


## 1. System requirements
### Hardware requirements
This code is requires only a standard computer with at least 8Gb of RAM to support typical single-cell RNA seq analyses.
### Software requirements
This code requires R version 4.2 or higher, and python code was written and tested with version 3.8.5. 

## 2. Installation guide
Installation should take less than 30 minutes.

### Python package installation
```{python}
conda create -n scvelo python=3.9
conda activate scvelo
pip install -U pandas anndata loompy scanpy scvelo igraph louvain pybind11 hnswlib
```

### R package installation
```{r}
# Install core package managers
install.packages(c(
  "renv",        # Reproducible environments
  "BiocManager", # Bioconductor installer
  "remotes",     # GitHub installs
  "devtools"     # Development utilities
))

# Install CRAN packages
install.packages(c(
  # Data structures and manipulation
  "Matrix", "dplyr", "data.table", "readr", "stringr", "reshape2",
  
  # Visualization and figure assembly
  "ggplot2", "cowplot", "patchwork", "ggrepel", "ggforce",
  "ggrastr", "pheatmap", "aplot", "ggdendroplot",
  "ggvenn", "ggVennDiagram",
  
  # Color palettes and scales
  "RColorBrewer", "paletteer", "MexBrewer", "scales",
  
  # Graphs and networks
  "igraph", "ggraph", "ggnetwork",
  
  # Statistical modeling
  "mgcv", "pcaPP",
  
  # Parallelization
  "future", "future.apply", "progressr",
  
  # Python interoperability
  "reticulate"
))

# Install Bioconductor packages
BiocManager::install(c(
  # Differential expression and normalization
  "edgeR", "DESeq2",
  
  # Single-cell infrastructure
  "Seurat", "SeuratDisk",
  
  # Signature and pathway scoring
  "SignatuR", "UCell", "AUCell", "escape", "fgsea",
  
  # Mixed models and tidy output
  "glmmTMB", "broom.mixed",
  
  # Trajectory and annotation
  "monocle3", "rtracklayer",
  
  # Parallel backends
  "BiocParallel"
))

# Install GitHub packages
# Note: install all dependent and suggested packages for seurat
remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", ref = "seurat5", quiet = TRUE)

# Note: RNA velocity analyses were performed using scVelo (Python), accessed via reticulate.
# Install the python packages first, then install the R package below
reticulate::use_condaenv("scvelo", required = TRUE)

```

## 3. Demo
There are no demo datasets provided for the code in this repository. Instead, the raw data can be downloaded and analyzed using the instructions in the following section. 

## 4. Instructions
Copy the GitHub repository.
```
git clone https://github.com/Prenauer/DKO-CAR-T
```
Download the processed data from GEO using the following links:
```
## Enter the following commands from the root directory of the GitHub repository.
wget -xP BULK_RNA/Data XXXX
wget -xP CRISPR_SCREEN/Data XXXX
```





