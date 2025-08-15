---
title: "sct_04_PathwayAnalysis"
author: "Paul Renauer"
date: "2025-08-14"
description: "This code performs generates (1) single-cell signatures of 
MSigDb Hallmark pathways, (2) runs differential pathway analyses of hallmark
signatures, (3) generates single-cell signatures of select phenotypes, and
(4) generates plots for all of these analyses."
output: html_document
---
```{r include = FALSE}
knitr::opts_chunk$set(message = F, eval=F)
```

### Load libraries and prepare environment
```{r Setup environment}
setwd('SCT/04_PathwayAnalysis')
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
library(escape)
library(ggraph)
library(igraph)
library(ggnetwork)
source('SCT/sct_00_functions.R')

## Set options
options(future.globals.maxSize= 10000*1024^2)
options(matrixStats.useNames.NA = 'deprecated')
options("ggrastr.default.dpi" = 750)

## Load SO
so <- readRDS('SCT/02_scVelo/Datasets/so_proc2.rds')

```

## Get single-cell signatures for pathways
```{r Get single-cell signatures for pathways}

## Prep for AUCell analysis
# Get Hallmark MSigDb pathways
gs <- getGeneSets(library='H')
# Get variable gene features
hvg <- VariableFeatures(so)


## Use escape package to run AUCell
es <- escape.matrix(so$RNA$counts[hvg,], gene.sets = gs, 
                    method='AUCell', groups = 10000, normalize=T,
                    min.size = 5,
                    BPPARAM = BiocParallel::MulticoreParam(workers = 20, 
                                                           force.GC=T, 
                                                           progressbar=T))
# Export AUCell signature matrix
saveRDS(es, 'escapeAUC_hmHVG.rds')


## Add hallmark pathway signatures as an assay object
so$es <- CreateAssay5Object(count=t(es), data=t(es))

```

### Differential pathway analysis
```{r Differential pathway analysis}

## Set assay and celltype groups
DefaultAssay(so) <- 'es'
groups <- c('TPEX/TEX','TMEM','TRM','Act.T','TE', 'CD4T')


## Differential pathway analysis (one-vs-all: determine highly-var. pathways)
# Iterate through groups
pw <- do.call(rbind, lapply(groups, function(celltype){
    # Run one-vs-all analysis
    de <- FindAllMarkers(subset(so, subset=(groups==celltype)), 
                         group.by='sample', logfc.threshold=0, densify=T)
    # Change column names
    colnames(de) <- c('pval','logFC','pct.1','pct.2', 'padj','group2','feature')
    # Add columns to table
    de <- data.frame(celltype=celltype, group1=de$group2, de,
                     row.names=NULL)
    de$group2 <- 'all'
    # re-order table
    de <- de[with(de, order(-log10(pval), abs(logFC), decreasing=T)),]
    # return
    return(de)
}))
# re-order table
pw <- pw[with(pw, order(celltype, group1, pval)),]
# Exclude pathways that are irrelevant to T cell tumor-immunity
pw <- pw[-grep('MYCOBACTERIUM-TUBERCULOSIS|HALLMARK-UV-RESPONSE|COAGULATION|SPERMATOGENESIS|MYOGENESIS', 
               pw$feature),]
# Export results
write.table(pw, 'pwHM_AucHvg_OneVsAll.txt', sep='\t', quote=F, row.names=F)


## Differential pathway analysis (pairwise)
genos <- c('DKO','NR4A1','SOCS3','GTC')
# Iterate through pairwise genotype comparisons
de <- do.call(rbind, lapply(1:3, function(i){
    do.call(rbind, lapply((i+1):4, function(j){
        gp1 <- genos[i]
        gp2 <- genos[j]
        # Iterate through groups
        do.call(rbind, lapply(groups, function(celltype){
            # Run DE analysis
            de <- FindMarkers(subset(so, subset=(groups==celltype)), 
                              logfc.threshold=0, densify=T,group.by='sample',
                              features=rownames(so$es), ident.1=gp1,ident.2=gp2)
            # Change column names
            colnames(de) <- c('pval','logFC','pct.1','pct.2', 'padj')
            # Add columns to table
            de <- data.frame(celltype=celltype, feature=rownames(de), 
                             comp=paste0(gp1,'-',gp2), group1=gp1, group2=gp2, de)
            # re-order table
            de <- de[with(de, order(-log10(pval), abs(logFC), decreasing=T)),]
            # return table
            return(de)
        }))
    }))
}))
# Export table
write.table(de, 'de_pathways_v1.0.txt', sep='\t', quote=F, row.names=F)

```

### Plot Volcanos of Diff. pathways
```{r Plot Volcanos of Diff. pathways}

## Load diff-pw results for variable pathways
de <- read.delim('pwHM_AucHvg_OneVsAll.txt')
# subset results for DKO-specific variable pathways
de <- de[which(de$group1=='DKO'),]


## Prep table for plots
# Add comparison column
de$comp <- 'DKO-vs-all'
comps <- unique(de$comp)
de$group <- de$comp # need the group column for the volcano plot function
# Exclude redundant part of pathway names
de$feature <- str_remove(de$feature, 'HALLMARK-')


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
ggsave(plot=p, 'Figures/vol_de_pathways_oneVsAll.pdf', height=1.25*6, 
       width=1.25*1, scale=2)


## Prep for pairwise volcano plots
# Load result table
de <- read.delim('de_pathways_v1.0.txt')
comps <- unique(de$comp)
de$group <- de$comp # need the group column for the volcano plot function
# Exclude redundant part of pathway names
de$feature <- str_remove(de$feature, 'HALLMARK-')

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
ggsave(plot=p, 'Figures/vol_de_pathways_pairwise.pdf', height=1.25*6, 
       width=1.25*6, scale=2)

```

### Make upset plots to compare Diff. pathway results
```{r Make upset plots to compare Diff. pathway results}

## Prep data for plots
# load pairwise results
de <- read.delim("de_pathways_v1.0.txt")
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
ggsave(plot=p, 'Figures/upset_de_pathways.pdf', height=1.5*2, 
       width=2*6, scale=2.2)

```

### Pathway network analyses
```{r Pathway network analyses}


## Prep data for plotting
# create id showing celltype and comparison
pw$group <- paste0(pw$celltype,': ',pw$group1,'-',pw$group2)
# Fix columns
pw$ct <- pw$celltype
pw$pathway <- str_replace_all(pw$feature,'-','_')
pw$NES <- pw$logFC
# Get separate input dataframes for up and down-regulated pathways
up <- slice_head(pw[which(pw$logFC > 0 & pw$padj < 0.01),], by='group',n=2)
dn <- slice_head(pw[which(pw$logFC < 0 & pw$padj < 0.01),], by='group',n=2)
# generate and combine network plots
p <- cowplot::plot_grid(RunNetPW(up, 'Upreg. pathways', 23), 
                        RunNetPW(dn, 'Downreg. pathways', 23),
                        align='hv', axis='lbt',ncol=2)
# Export plots
ggsave(plot=p, paste0('Figures/net_pwHM_AucHvg_OneVsAll_v0.8d.pdf'), 
       height=5*1, width=4*2, scale=1.3)


```

### Generate single-cell signatures for phenotypes
```{r Generate single-cell signatures for phenotypes}

## Create list of pathway gene ontologies from MSigDB (Reactome, Hallmark, & TF)
gmt <- list(fgsea::gmtPathways(
    'SCT/00_ReferenceInfo/c2.reactome_pathways.v7.4.symbols.gmt.txt'),
            fgsea::gmtPathways(
                'SCT/00_ReferenceInfo/MSigDB_Hallmark_2020.txt'),
    fgsea::gmtPathways(
                'SCT/00_ReferenceInfo/c3.TF_targets.gtrd.v7.4.symbols.gmt.txt')
    ) %>% unlist(.,F)
# Filter pathway genes by those with count data
gmt <- lapply(gmt, function(x) intersect(x, rownames(so)))


## Create gene set lists (both custom and from MSigDb)
pw.list <- list(list(OXPHOS=gmt[['Oxidative Phosphorylation']],
                     PGC1A_targets=gmt[['PPARGC1A_TARGET_GENES']]),
                list(IL2_STAT5=gmt[['IL-2/STAT5 Signaling']],
                     IFNa=gmt[['Interferon Alpha Response']],
                     TNF_NFkB=gmt[['TNF-alpha Signaling via NF-kB']],
                     mTORC1=gmt[['mTORC1 Signaling']],
                     MYC1=gmt[['Myc Targets V1']],
                     MYC2=gmt[['Myc Targets V2']]),
                list(Hypoxia=gmt[['Hypoxia']]),
                list(TRM=c('CD69','ZNF683','CXCR3','CXCR6'),
                     Exh=c('TIGIT','HAVCR2','TOX','ENTPD1'),
                     MEM=c('KLF2','SELL','EOMES','S1PR1','IL15RA'), #,
                     TPEX=c('SLAMF6','TOX'),
                     Cytotoxic=c('GZMA','GZMB','PRF1','FASLG'),
                     Cytokines=c('TNF','IFNG'),
                     Prolif=c('TOP2A','MKI67'))
) %>% unlist(.,F)


## Generate matrix of phenotype signatures using AUCell via escape package
e <- escape::escape.matrix(so$RNA$counts[VariableFeatures(so$RNA),], 
                           gene.sets = pw.list, 
                           method='AUCell', groups = 10000, #normalize=T,
                           min.size = 2,
                           BPPARAM = 
                               BiocParallel::MulticoreParam(workers = 20, 
                                                            force.GC=T, 
                                                            progressbar=T))
# export signature matrix
saveRDS(e, 'aucell_phenotypeScores_v1.0.rds')

```

### Violin plots to compare phenotype signatures
```{r Violin plots to compare phenotype signatures}


## Add signatures to SO
# Load phenotype signature data
e <- readRDS('aucell_phenotypeScores_v1.0.rds')[colnames(so),]
# Add scaled signature matrix to SO
so@meta.data <- data.frame(so@meta.data, scale(e))


## Prep data for plotting
# set colors
cols <- rep((c('#B2B2B2','#496AB4','#E9AC4C','#D85446')),6)
# order the groups and genotypes
so$groups <- factor(so$groups, levels=
                        c('TMEM','TRM','TE','Act.T','TPEX/TEX','CD4T'))
so$sample <- factor(so$sample, levels=c('GTC','NR4A1','SOCS3','DKO'))



## Make violin of signaling signatures
features <- c('IL2_STAT5','TNF_NFkB','mTORC1','OXPHOS', 'Hypoxia')
VlnPlot(so,features,group.by='groups', pt.size=0, stack=T, fill.by='ident', 
        adjust=0.6,split.by='sample',
        cols=cols, same.y.lims=T, flip=T) + ylim(c(-1,5)) +
    geom_hline(yintercept=1, linetype='dashed', color='gray10', linewidth=0.5)+
    labs(y='Gene signature (z score)', x=NULL) + 
    theme(legend.position='none',
          legend.title=element_text(size=10))
ggsave(plot=p, 'Figures/topPW_Vln_Signaling_v1.2.pdf', 
       height=3.2, width=6,scale=1.5)

## Make violin of Mito-Fitness signatures
features <- c('PGC1A_targets','PPARA_targets','mtBiogenesis')
p <- VlnPlot(so,features,group.by='groups', pt.size=0, stack=T, fill.by='ident', 
             adjust=0.6,split.by='sample',
             cols=cols, same.y.lims=T, flip=T) + ylim(c(0, 4)) + 
    geom_hline(yintercept=1, linetype='dashed', color='gray10', linewidth=0.5)+
    labs(y='Gene signature (z score)', x=NULL) + 
    theme(legend.position='none',
          legend.title=element_text(size=10))
ggsave(plot=p, 'Figures/topPW_Vln_Metabolism_v1.0.pdf', 
       height=2.5, width=6,scale=1.5)


## Make violin of phenotype signatures
features <- c('Prolif.1','Cytokines.1','Cytotoxic.1','Exh.1')
p <- VlnPlot(so,features,group.by='groups', pt.size=0, stack=T, fill.by='ident', 
             adjust=0.6,split.by='sample',
             cols=cols, same.y.lims=T, flip=T) + ylim(c(0,3)) +
    geom_hline(yintercept=1, linetype='dashed', color='gray10', linewidth=0.5)+
    labs(y='Gene signature (z score)', x=NULL) + 
    theme(legend.position='none', 
          legend.title=element_text(size=10))# +
ggsave(plot=p, 'Figures/topPW_Vln_Phenotype_v1.0.pdf', 
       height=2.5, width=6,scale=1.5)

```


### Get stats for violin plots
```{r Get stats for violin plots}


## Get stats for violin plots
# get data
d <- so@meta.data
# Iterate through pathways/phenotypes
anova_res <- do.call(rbind, lapply(
    c('IL2_STAT5','TNF_NFkB','mTORC1','OXPHOS', 'Hypoxia'), function(param){
        # subset by signature                                   
        d$y <- d[,param]
        # run 2-way anova
        anova <- aov(y ~ groups * sample, data = d)
        # run Tukey's post-hoc analysis
        r <- TukeyHSD(anova, which=c('groups:sample'))$`groups:sample` %>% 
            data.frame()
        # Reorganize result table
        tmp <- str_split(rownames(r), ':|-', simplify=T)
        colnames(tmp) <- c('c1','g1','c2','g2')
        r <- cbind(r,tmp)
        # Filter results to only intra- cellgroup comparisons
        r <- r[which(r$c1==r$c2),]
        # Create a column to order the results relative to plot 
        r$c1 <- factor(r$c1, levels=c('TMEM','TRM','TE','TPEX/TEX',
                                      'Act.T', 'CD4T'))
        r$g <- apply(r,1,function(x){
            paste0(sort(factor(x[c('g1','g2')], 
                               levels=c('GTC','NR4A1','SOCS3','DKO'))),
                   collapse='-')
        })
        r$g <- factor(r$g, levels=c('GTC-NR4A1','GTC-SOCS3','GTC-DKO',
                                    'NR4A1-SOCS3','NR4A1-DKO','SOCS3-DKO'))
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
write.table(anova_res, 'Data/vln_signature_stats_v1.0.txt', quote=F, 
            row.names=F, sep='\t')


```



