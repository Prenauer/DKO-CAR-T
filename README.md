# DKO CAR-T manuscript code


## 1. System requirements
### Hardware requirements
This code is requires only a standard computer with at least 8Gb of RAM to support typical single-cell RNA seq analyses.
### Software requirements
This code requires R version 4.2 or higher, and python code was written and tested with version 3.8.5. 

## 2. Installation guide
Installation should take less than 30 minutes.
### R package installation
```{r}
install.packages(c('renv','BiocManager','remotes','devtools','dplyr',
'stringr','patchwork','reshape2','factoextra','RColorBrewer','pheatmap',
'leidenbase'))
BiocManager::install(c('edgeR','zellkonverter','rhdf5','decoupleR','glmGamPoi'))
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
```
#Note: install all dependent and suggested packages for seurat

### Python package installation
```{python}
pip install -U pandas scvelo igraph louvain pybind11 hnswlib
```

## 3. Demo
There are no demo datasets provided for the code in this repository. Instead, the raw data can be downloaded and analyzed using the instructions in the following section. 

## 4. Instructions
Copy the GitHub repository.
```
git clone https://github.com/Prenauer/2024_DKOT
```
Download the processed data from GEO using the following links:
```
## Enter the following commands from the root directory of the GitHub repository.
wget -xP BULK_RNA/Data XXXX
wget -xP CRISPR_SCREEN/Data XXXX
```





