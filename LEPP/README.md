## LEPP: low gene expression as a predictor of T cell phenotype

Author: Paul Renauer

Date: 2026-01-14

Description: 
This uses public single-cell RNA-seq data, curated from the 
ProjecTILs tumor-infiltrating T lymphocytes database to assess whether the 
double-low expression of target genes is a better predictor of T cell function  
than single-low expression of these genes, using gene module scores to 
represent T cell phenotypes.
##

### Script contents

```{r }
###############################################################################
# LEPP analysis driver script
#
# This script:
#   1. Loads required packages
#   2. Sources LEPP function library
#   3. Loads and preprocesses data
#   4. Computes phenotype signatures
#   5. Runs LEPP single + pairwise models
#   6. Writes result tables to disk
#
# All statistical functions are found in: scripts/00_lepp_functions.R
# ProjecTILs data found at: https://figshare.com/articles/dataset/ProjecTILs_human_reference_atlas_of_CD8_tumor-infiltrating_T_cells_CD8_TIL_version_1/23608308?file=41414556
###############################################################################

```

### Load packages

```{r setup, message=FALSE, warning=FALSE}
# Load all required packages quietly
suppressPackageStartupMessages({

  # Sparse matrices and expression handling
  library(Matrix)              # Efficient sparse matrix support
  library(edgeR)               # Expression utilities

  # Data manipulation
  library(dplyr)               # Data wrangling
  library(readr)               # File I/O
  library(Seurat)              # Single-cell container

  # Mixed models
  library(glmmTMB)             # Gaussian mixed-effects models
  library(broom.mixed)         # Tidy extraction of mixed models

  # Signature scoring
  library(SignatuR)            # Curated gene signatures
  library(UCell)               # Rank-based signature scoring
  library(BiocParallel)        # Parallel backend for UCell

  # Parallel execution
  library(future)              # Parallel planning
  library(future.apply)        # Parallel lapply
  library(progressr)           # Progress handling

  # Plotting
  library(ggplot2)             # Grammar of graphics
  library(cowplot)             # Plot composition
  library(reshape2)            # Data reshaping
  library(scales)              # Color/scale helpers
  library(paletteer)           # Color palettes
  library(MexBrewer)           # Additional color palettes
})

# Increase allowable size for future globals (needed for large objects)
options(future.globals.maxSize = 10 * 1024^3)

# Silence deprecated matrixStats naming warnings
options(matrixStats.useNames.NA = "deprecated")

# Register multicore backend
plan(multicore, workers = 18)

# Disable nested progress bars
progressr::handlers(global = FALSE)

# Source LEPP statistical function library
source("scripts/00_lepp_functions.R")
```

## Set parameters

```{r params}

# Note: the seurat object can be downloaded here: https://figshare.com/articles/dataset/ProjecTILs_human_reference_atlas_of_CD8_tumor-infiltrating_T_cells_CD8_TIL_version_1/23608308?file=41414556
# Path to input Seurat object
seurat_rds_path <- "datasets/Projectils_Cd8T.rds"

# Output directories
results_dir <- "data"
figures_dir <- "figures"

# Create output directories if they do not exist
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Statistical thresholds
alpha_single <- 0.05        # FDR threshold for single-gene tests
alpha_pair   <- 0.05        # FDR threshold for pairwise tests

# Bootstrap parameters
n_boot      <- 2000         # Number of bootstrap samples
seed        <- 1            # Random seed
chunk_size <- 1000          # Bootstrap chunk size

# Phenotype signatures to analyze
phenotype_signatures <- c("Tcell.cytotoxicity","IFN","cellCycle.G2M",
  "Tcell.stemness","Tcell.exhaustion")
```

## Load and filter data

```{r data}
# Load Seurat object from disk
s <- readRDS(seurat_rds_path)

# Restrict to singlets with functional annotation
s <- subset(s,subset = (db.class == "singlet" & !is.na(functional.cluster)))
```

## Compute phenotype signatures

```{r signatures}
# Load bundled SignatuR data
data(SignatuR)

# Extract relevant human signatures
sig_list <- GetSignature(SignatuR$Hs)[phenotype_signatures]

# Compute UCell scores using parallel backend
ucell_scores <- ScoreSignatures_UCell(
  s$RNA$counts,
  features = sig_list,
  BPPARAM = MulticoreParam(workers = 18, progressbar = TRUE)
)

# Rename signature score columns for readability
colnames(ucell_scores) <- c("Cytotoxicity","IFN","Proliferation",
                            "Stemness","Exhaustion")

# Store phenotype names
phenos <- colnames(ucell_scores)
```

### Assemble analysis table

```{r table}
# Construct analysis data frame combining expression, phenotypes, and metadata
analysis_df <- data.frame(
  t(s$RNA$counts[pred, ]),                # Expression for predictor genes
  ucell_scores[colnames(s), ],            # Phenotype scores
  subset = factor(s$functional.cluster),  # Functional cluster
  patient_id = factor(s$Sample),          # Patient ID
  tissue = factor(s$Tissue),              # Tissue of origin
  nCount_RNA = s$nCount_RNA                # Library size
)
```

### Run LEPP analysis

```{r analysis}
# Set seed for reproducibility
set.seed(seed)

# Run LEPP scan across all phenotypes
lepp_results <- bind_rows(lapply(phenos, function(ph) {

  # Print progress message
  message("Running LEPP for phenotype: ", ph)

  # Attach phenotype as response variable
  df <- transform(
    analysis_df,
    phenotype_score = analysis_df[, ph]
  )

  # Run single + pairwise LEPP scan
  out <- run_scan_parallel(
    df = df, genes = pred, alpha_single = alpha_single, alpha_pair = alpha_pair,
    n_boot = n_boot, seed = seed, chunk_size = chunk_size)

  # Annotate phenotype
  transform(out, phenotype_score = ph)
}))
```

### Save analysis results

```{r save-results}
# Write LEPP results to disk
write.table(lepp_results, 
            file = file.path(results_dir, "lepp_results_v0.1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE )
```

### Figure: Gene dominance (pair vs best single)

```{r gene-dominance}
# Subset to pair-vs-single tests
d <- lepp_results[lepp_results$test_type == "pair_vs_single", ]

# Construct pair identifier
d$id <- apply(
  d[, c("geneA", "geneB")], 1, function(x) paste0(x[!is.na(x)], collapse = "_"))

# Select color palette
cols <- as.character(paletteer::paletteer_d("MetBrewer::Degas")[-1])

# Identify significant comparisons
d_sig <- d %>%
  filter(sig) %>%
  mutate(label = paste0("q=", formatC(q, format = "e", digits = 2)))

# Generate dot–whisker plot
p <- ggplot(d, aes(y = factor(id), x = beta, color = factor(phenotype_score))) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point() +
  geom_errorbar(
    aes(xmin = delta_ci_low, xmax = delta_ci_high),
    width = 0.5, size = 1) +
  geom_text(
    data = d_sig, aes(label = label), 
    position = position_nudge(y = 0.5), size = 3) +
  facet_wrap(vars(phenotype_score), scales = "free_x", ncol = 5) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(x = "Delta (pair − best single)", y = "Gene pair",
    color = "Phenotype",  title = "Gene dominance: pair-low vs best-single-low")

# Save figure
ggsave(plot = p, 
       filename = file.path(figures_dir, "lepp_gene_dominance_v0.1.pdf"),
       height = 3, width = 6, scale = 1.5)
```

### Figure: Single vs paired low-expression effects

```{r dotplot}
# Create identifier for single and paired tests
make_id <- function(a, b)
  apply(cbind(a, b), 1, function(x) paste0(x[!is.na(x)], collapse = "_"))

# Prepare plotting data
d <- lepp_results %>%
  filter(test_type != "pair_vs_single") %>%
  mutate(logq = -log10(q), id = make_id(geneA, geneB), pheno = phenotype_score)

# Restructure data for dot plot
d_plot <- do.call(rbind, lapply(unique(d$pheno), function(ph) {
  do.call(rbind, lapply(unique(d$id[d$test_type == "pair"]), function(id) {
    pair <- c(unlist(strsplit(id, "_")), NA)
    tmp <- d[d$geneA %in% pair & d$geneB %in% pair & d$pheno == ph, ]
    idx <- c(which(tmp$id == pair[1]), 
             which(tmp$id == pair[2]), which(tmp$id == id))
    labs <- c("A", "B", "AB")
    data.frame(pheno = ph, pair = id, x = labs[idx], tmp[, c("beta", "logq")])
  }))
}))

# Order x-axis
d_plot$x <- factor(d_plot$x, levels = c("A", "B", "AB"))

# Define color scale
col_fun <- colorRampPalette(
  rev(MexBrewer::MexPalettes[["Huida"]][[1]][2:10])
)

# Generate dot plot
p <- ggplot(d_plot, aes(x = x, y = pair)) +
  geom_point(aes(fill = beta, size = logq), shape = 21, stroke = 0.4) +
  scale_fill_gradientn(colors = col_fun(50)) +
  scale_size(range = c(2.5, 6)) +
  facet_wrap(vars(pheno), ncol = 5, scales = "free_x") +
  theme_classic() +
  labs( x = "Single vs paired low expression", y = "Gene pair",
    fill = "β", size = "-log10(q)", 
    title = "Effect of low gene expression on phenotypes")

# Save dot plot
ggsave(plot = p,
  filename = file.path(figures_dir, "lepp_dotplot_v0.1.pdf"),
  height = 2.5, width = 5, scale = 1.5 )
```
