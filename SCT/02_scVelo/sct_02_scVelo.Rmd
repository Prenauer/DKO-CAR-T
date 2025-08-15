---
title: "sct_02_scVelo"
author: "Paul Renauer"
date: "2025-08-14"
description: "This code includes steps to (1) prepare data for scVelo analysis
and (2) perform post-analyses for scVelo results."
output: html_document
---
```{r include = FALSE}
knitr::opts_chunk$set(message = F, eval=F)
```

### Load libraries and prepare environment
```{r Setup environment}
setwd('SCT/02_scVelo')
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrastr)
library(Matrix)
library(reticulate)
library(Seurat)
library(dplyr)
library(data.table)
library(SeuratDisk)
library(SeuratWrappers)
source('SCT/sct_00_functions.R')

## Set options
options(future.globals.maxSize= 10000*1024^2)
options(matrixStats.useNames.NA = 'deprecated')
options("ggrastr.default.dpi" = 750)
use_condaenv('py39_native')

## Load SO
so <- readRDS('SCT/01_Processing/Datasets/so_proc1.rds')

```

### Export data for scvelo
```{r Export data for scvelo}

## Rename cell ids to match those from velocyto outputs
# Map the genotype labels to the random string-ids from velocyto
id <- c(DKO='U7O49', GTC='TDCRY', NR4A1='KT7Z1', SOCS3='01U4K')
# Create new cell ids
new.names=paste0('possorted_genome_bam_',id[so$sample],':',
                 str_split_i(colnames(so), '-', 1),'x')
# Rename cells
so2 <- RenameCells(so, new.names=new.names)


## Subset data to CD8 T cells only, so trajectories are not incorrectly 
##    drawn between CD4 and CD8 T cells.
so2 <- subset(so2, subset=(groups != 'CD4T'))


## Export loom data and the dimensional reduction embeddings
SaveLoom(so2, 'Datasets/so_proc2_namechange.loom', overwrite=T)
write.table(so2$umap@cell.embeddings, 
            'Datasets/umap_embedding.dat', sep='\t',
            quote=F, row.names=F)
write.table(so2$inmfNorm@cell.embeddings, 
            'Datasets/nmf_embedding.dat', sep='\t',
            quote=F, row.names=F)
rm(so2)

```

### Import scVelo analysis results
```{r Import scVelo analysis results}

## Import scVelo observation data
velo.meta <- read.delim('Data/04_scvelo_obs_v2.3.txt')


## Convert cell ids back to that of SO dataset
velo.meta$id <- str_split_i(velo.meta$CellID, ':',2) %>% str_remove(., 'x') %>% 
    paste0(.,'-', structure(c(1:4), names=c('DKO','GTC','NR4A1','SOCS3'))[velo.meta$sample])
rownames(velo.meta) <- velo.meta$id


## Add scVelo result data to SO metadata
# Rename scVelo result columns
velo.meta <- velo.meta[colnames(so), 
                       c('velocity_pseudotime','velocity_self_transition',
                         'velocity_confidence_transition', 'latent_time', 
                         'velocity_length', 'velocity_confidence')]
# Add scVelo result data to SO metadata
so@meta.data <- cbind(so@meta.data, velo.meta)


## Export SO dataset
saveRDS(so, 'Datasets/so_proc2.rds')

```

### Plot results of scVelo analysis
```{r Plot results of scVelo analysis}

## Plot results of scVelo analysis 
# Set plotting colors
cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                  names=c('GTC','SOCS3','NR4A1','DKO'))
# Set parameters to plot
params <- c('velocity_pseudotime','velocity_self_transition',
            'latent_time','velocity_length','Prolif.1')
# Create list of beeswarm plots for each parameter
plist <- lapply(params, function(param){
    d <- so@meta.data[which(so$groups !='TH1'),]
    d$y <- d[,param]
    ggplot(d, aes(x=factor(groups, 
                           levels=rev(c('TMEM','TRM','TE','TPEX','Act.T'))), 
                         y=y, group=(sample))) +
    (geom_quasirandom_rast(dodge.width=0.8, color='gray20', size=0.8)) +
    (geom_quasirandom_rast(dodge.width=0.8, aes(color=(sample)), size=0.05)) +
    scale_color_manual(values=cols) + coord_flip() +
    labs(y=param, x='Cell type') +
    theme_classic()
})
# Combine plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv', axis='lbt', nrow=1)
ggsave(plot=p, 'Figures/bswarm_veloParams_v1.3.pdf', 
                 height=4,width=2.75*length(params))

```

## Perform stats for scVelo plots above
```{r Perform stats for scVelo plots above}

## For each scVelo parameter, perform 2-way anova with Tukey post-hoc analysis
anova_res <- do.call(rbind, lapply(params, function(param){
    # Create data frame with parameter data, making sure to omit the CD4T 
    d <- so@meta.data[which(so$groups !='CD4T'),]
    d$y <- d[,param]
    # Perform two-way ANOVA
    anova <- aov(y ~ groups * sample, data = d)
    # Perform Tukey post-hoc analyses
    r <- TukeyHSD(anova, which=c('groups:sample'))$`groups:sample` %>% 
        data.frame()
    # Reorganize result table
    tmp <- str_split(rownames(r), ':|-', simplify=T)
    colnames(tmp) <- c('c1','g1','c2','g2')
    r <- cbind(r,tmp)
    # Filter results to only intra- cellgroup comparisons
    r <- r[which(r$c1==r$c2),]
    # Create a column to order the results relative to plot 
    r$c1 <- factor(r$c1, levels=c('TMEM','TRM','TE','TPEX/TEX','Act.T'))
    r$g <- apply(r,1,function(x){
        paste0(sort(factor(x[c('g1','g2')], levels=c('GTC','NR4A1','SOCS3','DKO'))),collapse='-')
    })
    r$g <- factor(r$g, levels=c('GTC-NR4A1','GTC-SOCS3','GTC-DKO','NR4A1-SOCS3','NR4A1-DKO','SOCS3-DKO'))
    r <- r[with(r, order(c1,g)),]
    # Create column with astrices to represent significance
    r$sig <- sapply(r$p.adj, function(x) {
        case_when(x < 1e-4 ~ '****',
                  x < 1e-3 ~ '***',
                  x < 1e-2 ~ '**',
                  x < 0.05 ~ '*',
                  x >= 0.05 ~ 'ns')
    })
    # return results
    return(cbind(variable=param, r))
}))
# Export table
write.table(anova_res, 'Data/dk_comp_groups-geno_v1.0.txt', quote=F, row.names=F, sep='\t')

```

### Make plots for top diff. kinetic genes
```{r Make plots for top diff. kinetic genes}

## Select top kinetic genes
# Load differential kinetics (DK) results 
d <- read.delim('Data/04_scvelo_var_v2.3.txt')
# Adjust p values for multiple testing
d$fit_padj_kinetics <- p.adjust(d$fit_pval_kinetics, n=nrow(d))
# Filter genes for significance and quality model fitting
d <- d[which(d$fit_padj_kinetics < 0.01 & d$spearmans_score > 0.5),]
# make vector of top genes
dkg <- d$Gene


## Make heatmap and PCA of top DK genes
p <- cowplot::plot_grid(plotlist=MakeDotPlot(dkg, 'DK Genes', 1/3,'ward.D2'), 
                        align='h', axis='lbt',ncol=2, rel_widths=c(0.5,0.5))
ggsave(plot=p, 'velo_dkgenes_hm_pca_v1.0.pdf', limitsize = F, height=4,width=6.5, scale=1.5)

```






