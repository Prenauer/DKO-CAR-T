---
title: "sct_03_DE-Genes"
author: "Paul Renauer"
date: "2025-08-14"
description: "This code performs differential expression analyses of genes
and generates plots."
output: html_document
---

```{r include = FALSE}
knitr::opts_chunk$set(message = F, eval=F)
```

### Load libraries and prepare environment

```{r Setup environment}
setwd('SCT/03_DE_Genes')
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
library(ggvenn)
library(ggVennDiagram)
source('SCT/sct_00_functions.R')

## Set options
options(future.globals.maxSize= 10000*1024^2)
options(matrixStats.useNames.NA = 'deprecated')
options("ggrastr.default.dpi" = 750)

## Load SO
so <- readRDS('SCT/02_scVelo/Datasets/so_proc2.rds')

```

### Run pairwise DE analyses

```{r Run pairwise DE analyses }

## Prep for DE analyses
# Set list of genes with include in DE analyses
genes2use <- grep('^HIST|^AC([0-9])|^AL([0-9])|^AP([0-9])|^AF([0-9])|^CU([0-9])|^RPL([0-9])|^RPS([0-9])|^LINC|orf', 
                  rownames(so), invert=T, value=T)
# set celltype-groups and genotypes to interate through.
groups <- c('TPEX','TMEM','TRM','Tcell','TE','TH1')
genos <- c('DKO','NR4A1','SOCS3','GTC')


## Run pairwise DE analyses
de <- do.call(rbind, lapply(1:3, function(i){
    do.call(rbind, lapply((i+1):4, function(j){
        # set group names
        gp1 <- genos[i]
        gp2 <- genos[j]
        # Iterate through celltype-groups
        do.call(rbind, lapply(groups, function(celltype){
            # DE analysis
            de <- FindMarkers(subset(so, subset=(groups==celltype)), 
                              logfc.threshold=0, densify=T, group.by='sample',
                              features=genes2use, ident.1=gp1, ident.2=gp2)
            # Rename result columns
            colnames(de) <- c('pval','logFC','pct.1','pct.2', 'padj')
            # Create an output dataframe
            de <- data.frame(celltype=celltype, feature=rownames(de), 
                             comp=paste0(gp1,'-',gp2), group1=gp1, 
                             group2=gp2, de)
            # Re-order the result table
            de <- de[with(de, order(-log10(pval), abs(logFC), decreasing=T)),]
            # Return dataframe
            return(de)
        }))
    }))
}))

## Export table
write.table(de, 'Data/de_genes.txt', sep='\t', quote=F, row.names=F)

```

### Make volcano plots for pairwise DE analyses

```{r Make volcano plots for pairwise DE analyses}

# Iterate through genotype comparisons
plist <- lapply(1:3, function(i){
    lapply((i+1):4, function(j){
        # set group names
        gp1 <- genos[i]
        gp2 <- genos[j]
        # Iterate through celltype-groups
        lapply(groups, function(celltype){
            # Make a volcano plot
            plotVolcano(de[which(de1$group1==gp1 & de1$group2=='GTC' & 
                                     de1$celltype==celltype),], 
                        pt.size=0.8, rastr=T, lfc.thresh=1.5, 
                        title=paste0(celltype,': ',gp1,'-GTC'))
        })
    }) %>% unlist(recursive=F)
}) %>% unlist(recursive=F)
        
## Combine to one plot and export       
p <- cowplot::plot_grid(plotlist=plist, align='hv',axis='lbt',ncol=6, byrow=F)
ggsave(plot=p, 'Figures/vol_de_genes.pdf', height=2.5*6, width=2.5*6)

```

### Make upset plots to compare DE genes across celltypes and genotypes

```{r Make upset plots to compare DE genes across celltypes and genotypes}

## Create List of intersections between DEGs of each comparison
# Label up and down regulation
de$dir <- c('up','dn')[1+as.integer(de$logFC < 0)]
# Label significance
de$sig <- (de$padj < 0.01 & abs(de$logFC) > 1.5)
# Generate unique ids for DEG categories
de$id <- paste0(de$dir, ':', de$celltype, ':', de$comp)
# get gene list
g <- unique(de$feature)
# Generate 0/1 table describing which genes were signif in which DEG category
d <- do.call(cbind, lapply(unique(de$id), function(id){
    g %in% de$feature[which(de$id==id & de$sig)]
}))
colnames(d) <- unique(de$id)
# Create list of intersection info
d <- ggvenn::data_frame_to_list(data.frame(d))
# Fix names of list
names(d) <- str_replace(names(d), 'up.', 'UP_') %>% str_replace(., 'dn.', 'DN_') %>% 
    str_replace(., '\\.', ':') %>% str_replace(., '\\.', '-')


## Upset plots
# Iterate through celltype groups
plist <- lapply(groups, function(ct){
    ## make venn object
    v <- Venn(d[grep(ct,names(d))])
    
    # get total number of intersections > 0
    n.intersects <- sum(process_region_data(v)$count > 0)
    
    ## get colors
    # get info from venn object
    rd <- process_region_data(v)
    rd <- str_split(head(rd$name[order(-rd$count, rd$id)],n.intersects),
                    '/') %>% unlist()
    cols <- structure(c('firebrick','#496AB4'),names=c('UP','DN'))
    dot.col <- cols[str_split_i(rd, '_', 1)]
    bar.col <- cols[str_split_i(v@names, '_', 1)]
    
    
    ## Make upset plots
    # upset plots
    p <- plot_upset(v, n.intersects, sets.bar.color=bar.col, 
                    order.set.by='name',
                    relative_height=0.5,sets.bar.x.label='# Signif. genes',
                    intersection.matrix.color=dot.col) 
    # Convert plot to grob and return
    p <- cowplot::ggdraw() + cowplot::draw_grob(aplot::gglistGrob(p))
    return(p)
})
## Combine upset plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv',
                        axis='lbt',ncol=3)
ggsave(plot=p, 'Figures/upset_de_genes.pdf', height=1.5*2, width=2*6, scale=2.2)

```
