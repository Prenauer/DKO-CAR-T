---
title: "sct_05_TF_Analysis"
author: "Paul Renauer"
date: "2025-08-14"
description: "This code performs generates (1) single-cell signatures of 
MSigDb transcription facor binding sites, (2) runs differential signature 
analyses of the TF signatures, and (3) generates plots for all of these 
analyses."
output: html_document
---
```{r include = FALSE}
knitr::opts_chunk$set(message = F, eval=F)
```

### Load libraries and prepare environment
```{r Setup environment}
setwd('SCT/05_TF_Analysis')
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrastr)
library(Matrix)
library(reticulate)
library(Seurat)
library(dplyr)
library(data.table)
library(ggvenn)
library(ggVennDiagram)
library(AUCell)
library(aplot)
library(ggdendroplot)
library(BiocParallel)
library(mgcv)


source('SCT/sct_00_functions.R')

## Set options
options(future.globals.maxSize= 10000*1024^2)
options(matrixStats.useNames.NA = 'deprecated')
options("ggrastr.default.dpi" = 750)

## Load SO
so <- readRDS('SCT/02_scVelo/Datasets/so_proc2.rds')

```

## Get single-cell signatures for TFBS
```{r Get single-cell signatures for TFBS}

## Create TF signature scores
# Get highly variable genes
hvg <- VariableFeatures(so$RNA)
# build ranking dataset for AUCell
m <- AUCell_buildRankings(
    exprMat=so$RNA$counts[hvg,], featureType = "genes", splitByBlocks = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(tasks=100,force.GC=T,progressbar=T),
    verbose = T)
# Load gene set list for TFBSs
gs <- fgsea::gmtPathways(
                'SCT/00_ReferenceInfo/c3.TF_targets.gtrd.v7.4.symbols.gmt.txt')
# Assess AUCell signature matrix
m <- AUCell_calcAUC(gs, m, nCores=20, verbose=T)
m <- m@assays@data$AUC
# Export signature matrix
saveRDS(m, 'escapeAUC_tfHVG.rds')
rm(m)

```

### Differential TF-signature analysis 
```{r Differential TF-signature analysis}

## Add TF signatures to SO
# Load data
es.tf <- readRDS('escapeAUC_tfHVG.rds')
# Add TF signatures to SO as a new assay
so$es.tf <- CreateAssay5Object(count=t(es.tf), data=t(es.tf))
DefaultAssay(so) <- 'es.tf'


## Differential TF-signature between DKO and GTC of each celltype group
# Make vector of groups
groups <- c('TPEX/TEX','TMEM','TRM','Act.T','TE', 'CD4T')
# Iterate through celltype groups
de <- do.call(rbind, lapply(groups, function(celltype){
    # Run DKO-vs-GTC DE for a given celltype group
    de <- FindMarkers(subset(so, subset=(groups==celltype)), group.by='sample',
                      return.thresh=1, min.pct=0.001, ident.1='DKO', 
                      ident.2='GTC', logfc.threshold=0, densify=T)
    # Add columns to result table
    de <- data.frame(celltype=celltype de, geno='DKO',gene=rownames(de))
    # Rename columns of the result table
    colnames(de) <- c('celltype','pval','logFC','pct.1','pct.2', 
                      'padj', 'geno','gene')
    # Re-order table and return
    de <- de[with(de, order(-log10(pval), abs(logFC), decreasing=T)),]
    return(de)
    #}))
}))
# Export table
write.table(de,'Data/de_tf_results_v1.6.txt', sep='\t', quote=F, row.names=F)


## Pairwise differential TF-signature analysis between genotypes for groups.
genos <- c('DKO','NR4A1','SOCS3','GTC')
# Iterate through genotype comparisons
de <- do.call(rbind, lapply(1:3, function(i){
    do.call(rbind, lapply((i+1):4, function(j){
        gp1 <- genos[i]
        gp2 <- genos[j]
        # Iterate through celltype groups
        do.call(rbind, lapply(groups, function(celltype){
            # Subset SO by group then run DE
            de <- FindMarkers(subset(so, subset=(groups==celltype)),
                               group.by='sample',logfc.threshold=0, densify=T,
                              ident.1=gp1, ident.2=gp2,
                              features=rownames(so$es.tf))
            # Adjust result table and return
            colnames(de) <- c('pval','logFC','pct.1','pct.2', 'padj')
            de <- data.frame(celltype=celltype, feature=rownames(de), 
                             comp=paste0(gp1,'-',gp2), group1=gp1, group2=gp2, de)
            de <- de[with(de, order(-log10(pval), abs(logFC), decreasing=T)),]
            return(de)
        }))
    }))
}))
# Export table 
write.table(de, 'Data/de_tf_pairwise_v1.0.txt', sep='\t', quote=F, row.names=F)

```


### Plot Volcanos of Diff. TFs
```{r Plot Volcanos of Diff. TFs}

## Load diff-pw results for variable pathways
de <- read.delim('de_tf_results_v1.6.txt')
# subset results for DKO-specific variable pathways
de <- de[which(de$geno=='DKO'),]


## Prep table for plots
# Add comparison column
de$comp <- 'DKO-vs-all'
comps <- unique(de$comp)
de$group <- de$comp # need the group column for the volcano plot function
# Exclude redundant part of signature names
de$feature <- str_replace(de$feature, '-TARGET-GENES','_TF')


## Make volcano plots
# Iterate through groups
plist <- lapply(groups, function(celltype){
    # Iterate through comparisons
    lapply(comps, function(comp){
        # Make a plot title
       plot.title <- paste0(celltype,':DKO-vs-all')
        # Make volcano plot
        plotVolcano(de[which(de$celltype==celltype),], 
                    text.size=2, xlim.expansion=1.5,
                    pt.size=1, rastr=T, lfc.thresh=0.25, title=plot.title)
    })
}) %>% unlist(recursive=F)
# Combine plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv',
                        axis='lbt',ncol=length(comps))
ggsave(plot=p, 'Figures/vol_de_tf_oneVsAll.pdf', height=1.25*6, width=1.25*1, scale=2)


## Prep for pairwise volcano plots
# Load result table
de <- read.delim('de_tf_pairwise_v1.0.txt')
comps <- unique(de$comp)
de$group <- de$comp # need the group column for the volcano plot function
# Exclude redundant part of signature names
de$feature <- str_replace(de$feature, '-TARGET-GENES','_TF')

## Make volcano plots
# Iterate through groups
plist <- lapply(groups, function(celltype){
    # Iterate through comparisons
    lapply(comps, function(comp){
        # Make a plot title
        plot.title <- paste0(celltype,':',comp)
        # Make volcano plot
        plotVolcano(de[which(de$comp==comp & de$celltype==celltype),], 
                    text.size=2, xlim.expansion=3,
                    pt.size=1, rastr=T, lfc.thresh=0.25, title=plot.title)
    })
}) %>% unlist(recursive=F)
# Combine plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv',
                        axis='lbt',ncol=length(comps))
ggsave(plot=p, 'Figures/vol_de_tf_pairwise.pdf', height=1.25*6, 
       width=1.25*6, scale=2)

```

### Make upset plots to compare Diff. TF signature results
```{r Make upset plots to compare Diff. pathway results}

## Prep data for plots
# load pairwise results
de <- read.delim("Data/de_tf_pairwise_v1.0.txt")
# make columns for significance, direction, and id
de$dir <- c('up','dn')[1+as.integer(de$logFC < 0)]
de$sig <- (de$padj < 0.01 & abs(de$logFC) > 0.25)
de$id <- paste0(de$dir, ':', de$celltype, ':', de$comp)
# Get gene list
g <- unique(de$feature)
# Generate 0/1 table describing which genes were signif in which DEG category
d <- do.call(cbind, lapply(unique(de$id), function(id){
    g %in% de$feature[which(de$id==id & de$sig)]
}))
colnames(d) <- unique(de$id)
# Create list of intersection info
d <- ggvenn::data_frame_to_list(data.frame(d))
# Fix names of list
names(d) <- str_replace(names(d), 'up.', 'UP_') %>% 
    str_replace(., 'dn.', 'DN_') %>% 
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
                    relative_height=0.5,sets.bar.x.label='# Signif. pathways',
                    intersection.matrix.color=dot.col) 
    # Convert plot to grob and return
    p <- cowplot::ggdraw() + cowplot::draw_grob(aplot::gglistGrob(p))
    return(p)
})
# Combine plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv',
                        axis='lbt',ncol=3)
ggsave(plot=p, 'Figures/upset_de_tf.pdf', height=1.5*2, 
       width=2*6, scale=2.2)

```


### Make Dot plot heat maps for top results
```{r Make Dot plot heat maps for top results}

## Make plots for genes of top signatures
DefaultAssay(so) <- 'RNA'
plist <- list(MakeDotPlot(gs[['PPARGC1A_TARGET_GENES']], 'PPARGC1A TF Targets', 0.25,'ward.D2'),
              MakeDotPlot(gs[['PPARA_TARGET_GENES']], 'PPARA TF Targets', 0.25,'ward.D2'),
              MakeDotPlot(gs[['EPC1_TARGET_GENES']], 'EPC1 TF Targets', 0.25,'ward.D2'),
              MakeDotPlot(gs[['FOXR2_TARGET_GENES']], 'FOXR2 TF Targets', 0.25,'ward.D2')
              )
## Combine plots and export
p <- cowplot::plot_grid(plotlist=unlist(plist,F), align='h', axis='lbt',ncol=3, rel_widths=c(0.6,0.05,0.4))
ggsave(plot=p, 'Figures/hm_pca_TFsig_Genes_v1.0.pdf', limitsize = F, height=2.25*4,width=5.8, scale=3)

```

## Run Dynamic Signature Relationship (DSR) analyses between TFs and phenotypes
```{r Run Dynamic Signature Relationship (DSR) analyses between TFs and phenotypes}

## Setup TF data to run DSR
# Load diff-TF result table
de <- read.delim('de_tf_results_v1.5.txt')
# Get Signif. diff-TFs
detf <- de[which(de$geno=='DKO' & abs(de$logFC) > 0.25 & de$padj < 0.01), 
           c('celltype','gene')] %>% unique()
# Get matrix of TF signatures
m <- t(so$es.tf$data)[, unique(detf$gene)]


## Setup phenotype data to run DSR
# List phenotypes of interest
phenos <- c('MYC1','MYC2','IL2_STAT5','PGC1A_targets','OXPHOS','TNF_NFkB',
            'mTORC1','Hypoxia','Cytotoxic','Cytokines','Prolif', 'Exh')
# Load phenotype signature data
e <- readRDS('aucell_phenotypeScores_v1.0.rds')


## Create dataframe of batches to run for DSR
blist <- do.call(rbind, lapply(unique(detf$gene), function(g){
    do.call(rbind, lapply(phenos, function(pheno){
        data.frame(pheno=pheno, g=g)
    }))
}))


## Run batch GAMs
# Create blocks to run the DSR analyses
block_size <- 1000
blocks <- lapply(1:ceiling(nrow(blist)/block_size), function(i) 
    (((i-1)*block_size)+1):min((i*block_size), nrow(blist)))
# make dataframe to hold results
res <- data.frame()
# iterate blocks
for(j in 1:length(blocks)){
    block <- blocks[[j]]
    # Iterate through each block
    tmp <- bplapply(block, function(i){
        # Make input data frame
        d <- data.frame(tf=m[, blist$g[i]], pheno=e[, blist$pheno[i]], 
                        ct=factor(so$ct), 
                        s=factor(so$sample,
                                 levels=c('GTC','NR4A1','SOCS3','DKO')))
        # Filter input data
        d <- d[which(d$pheno > 0 & d$tf > 0),]
        # Run generalized additive model
        fit <- try(mgcv::bam(pheno ~ ct + s + s(tf, bs="ps",k=10),  
                              control=list(keepData=F),
                              data=d, Gamma(link='log'), method='fREML'), 
                    silent=T)
        if(attempt::is_try_error(fit)) return(NULL)
        # Get terms of model
        pr <- (predict(fit, newdata=d, type='terms'))
        # Get partial residuals for independent TF effect
        d$p_g <- (pr[,3] + residuals(fit))[rownames(d)]
        # Get partial residuals for TF+Genotype effect
        d$p_sg <- (rowSums(pr[,-2]) + residuals(fit))[rownames(d)]
        # run correlation test between TF and the TF effect on phenotype
        ctg <- cor.test(d$tf,d$p_g)[c('p.value','estimate')]
        # run correlation test between TF and the TF+Geno effect on phenotype
        ctsg <- cor.test(d$tf,d$p_sg)[c('p.value','estimate')]
        # Compile results into a dataframe
        ct <- structure(c(ctg, ctsg),
                        names=c('cor_g_p','cor_g_r2','cor_sg_p','cor_sg_r2'))
        r <- structure(c(summary(fit)$s.table[1,3:4],summary(fit)$r.sq, 
                         summary(fit)$dev.expl) , 
                       names=c('fit_F','fit_p','fit_r2','fit_dev'))
        res <- matrix(unlist(c(unlist(blist[i,]), ct, r)), nrow=1) %>% 
            data.frame()
        colnames(res) <- c(colnames(blist), names(ct), names(r))
        res[,3:ncol(res)] <- apply(res[,3:ncol(res)], 2, as.numeric)
        # Return results
        return(res)
    }, BPPARAM=MulticoreParam(workers=12, force.GC=F, progressbar=T,RNGseed=1))
    # Remove null results
    tmp <- tmp[!unlist(lapply(tmp,is.null))]
    # combine results into dataframe
    tmp <- do.call(rbind, tmp)
    # add block results to full results
    res <- rbind(res, tmp)
    # Collect garbage
    gc()
}
## Calculate adjusted p values
res <- mutate(res, .by=c('pheno'), 
              cor_g_padj=p.adjust(cor_g_p, n=length(rownames(so$RNA))),
              cor_sg_padj=p.adjust(cor_sg_p, n=length(rownames(so$RNA))),
              mod_padj=p.adjust(mod_p, n=length(rownames(so$RNA))))
## sort table and export
res <- res[order(res$pheno, -res$cor_g_r2),]
write.table(res, 'Data/DSR_tf-pheno_v1.8.txt', sep='\t', quote=F, row.names=F)

```

### Heatmap of DSR results
```{r Heatmap of DSR results}

## list phenotypes of interest
phenos <- c('MYC1','MYC2','IL2_STAT5','OXPHOS','TNF_NFkB',
            'mTORC1','TNFR2_ncNFkB','Hypoxia','Prolif')
# Load DSR results
res <- read.delim('Data/DSR_tf-pheno_v1.8.txt')
# filter results
targets <- unique(res$tf[which(res$mod_padj < 0.01 & 
                                  res$mod_r2 > 0.2 & 
                                   res$pheno %in% phenos)])
r <- res[which(res$tf %in% targets & res$pheno %in% phenos),]

# prep for plot
r <- reshape2::dcast(r, tf~pheno, value.var='cor_g_r2')
r <- data.frame(r[,-1], row.names=r[,1])
# Reformat colors so that it is winsorized to the minimum absolute extreme
library(RColorBrewer)
cols <- data.frame(cols = colorRampPalette(rev(
    brewer.pal(n = 6, name = "RdYlBu")))(100),
                   n=seq(-max(abs(r)), max(abs(r)), 2*max(abs(r))/99))
cols <- cols$cols[which(cols$n >= min(r) & cols$n <= max(r))]

# Make plot title
hm.title <- paste0('NR4A1_target-pheno map')
# Fix names
rownames(r) <- str_replace(rownames(r),'-TARGET-GENES', '_TF')
# Make plot and export
p <- cowplot::ggdraw() + cowplot::draw_grob(
    pheatmap::pheatmap(r,clustering_method='ward.D2', 
                       angle_col='90', fontsize=8,fontsize_row=8,treeheight_row=10, 
                       treeheight_col=10, cutree_cols=3, cutree_rows=2, 
                       main=hm.title)$gtable)
ggsave(plot=p, 'Figures/hm_DSR_tf-pheno_v1.8.pdf',  
       width=2.5, height=2, scale=1.5)

```

### Plot DSR results: TF signature vs TF-effects on Phenotype signature
```{r Plot DSR results: TF signature vs TF-effects on Phenotype signature}

## Create a dataframe where each line has the input of a plot 
r <- do.call(rbind, lapply(c('mTORC1','MYC1','OXPHOS'), function(pheno){
    do.call(rbind, lapply(paste0(c('PPARGC1A','PPARA'), '-TARGET-GENES'), 
                          function(g){
        do.call(rbind, lapply(c('TMEM','TRM','TE','TPEX/TEX','CD4T'), 
                              function(ct){
            data.frame(g=g, ct=ct, pheno=pheno)
        }))
    }))
}))
## Make scatter plots of TF signatures vs TF-effects on phenotype
plotlist=BiocParallel::bplapply(1:nrow(r), function(i){
    pheno <- r$pheno[i]
    g=r$g[i]
    celltype=r$ct[i]
    title <- paste0(celltype, ': ', pheno, '~',g)
    GamPlot(g, pheno, title, celltypes=celltype)
},BPPARAM = BiocParallel::MulticoreParam(workers = 10,  progressbar=T))
## Combine plots and export
p <- cowplot::plot_grid(plotlist=plotlist, align='hv', axis='lbt', ncol=10)
ggsave(plot=p, 'Figures/gamplots_DSR_tf-pheno_splitbyCelltype_v1.8.pdf', 
       height=2*3, width=2*10,scale=1.2)


## Make scatter plots of PPARGC1A TF signature vs adjusted effect on pathways
plist <- list(GamPlotTF(x='PPARGC1A-TARGET-GENES','OXPHOS'),
              GamPlotTF(x='PPARGC1A-TARGET-GENES','mTORC1'),
              GamPlotTF(x='PPARGC1A-TARGET-GENES','MYC1'),
              GamPlotTF(x='PPARGC1A-TARGET-GENES','MYC2'),
              GamEffectPlotTF(x='PPARGC1A-TARGET-GENES','OXPHOS'),
              GamEffectPlotTF(x='PPARGC1A-TARGET-GENES','mTORC1'),
              GamEffectPlotTF(x='PPARGC1A-TARGET-GENES','MYC1'),
              GamEffectPlotTF(x='PPARGC1A-TARGET-GENES','MYC2')
)
## Combine plots and export
p <- cowplot::plot_grid(plotlist=plist, align='hv', axis='lbt', ncol=8)
ggsave(plot=p, 'Figures/gamplots_DSR_tf-pheno_PGC1A.pdf', height=2*1, width=2*8,scale=1.2)

```
