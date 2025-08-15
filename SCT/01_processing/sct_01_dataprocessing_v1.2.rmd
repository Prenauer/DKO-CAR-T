---
title: "sct_01_Processing"
author: "Paul Renauer"
date: "2025-08-14"
description: "This code includes data preprocessing and processing steps.
Briefly, datasets are filtered, integrated, and reduced to 2-dimensions by UMAP.
Cells are then annotated and population sizes are compared by summary-level 
statistics."
output: html_document
---
```{r include = FALSE}
knitr::opts_chunk$set(message = F, eval=F)
```

### Load libraries and prepare environment
```{r Setup environment}
setwd('SCT/01_Processing')
library(dplyr)
library(stringr)
library(ggplot2)
library(Matrix)
library(reticulate)
library(Seurat)
library(dplyr)
library(data.table)
library(rliger)
source('SCT/sct_00_functions.R')

## Set options
options(future.globals.maxSize= 10000*1024^2)
options(matrixStats.useNames.NA = 'deprecated')
options("ggrastr.default.dpi" = 750)
use_condaenv('py39_native')

```

### Prepare Seurat Object
```{r Prepare Seurat Object}

# Set paths for raw data
rawDataPaths <- dir('cellranger', recursive=F, full.names=T)
rawDataPaths <- paste0(rawDataPaths, '/filtered_feature_bc_matrix')

# Read and merge matrices
m <- Read10X(rawDataPaths[1])
for(i in 2:4){
    rawDataPath <- rawDataPaths[i]
    m <- combineMatrices(m, Read10X(rawDataPath), paste0('-',i))
}


# Create Seurat object (so)
so <- CreateSeuratObject(counts=m, min.cells = 3, min.features=200)
rm(m)

# Add info to so
so$Percent.mt <- PercentageFeatureSet(so, pattern = '^MT-')
so$Low_Quality <- PercentageFeatureSet(so, features = c('MALAT1', 'KCNQ1OT1'))
groups <- c('DKO', 'GTC','NR4A1','SOCS3')
so$sample <- groups[as.integer(str_split_i(colnames(so), '-', 2))]

# Filter cells by mt % and quality
so <- subset(so, subset = Percent.mt < 20 & Low_Quality < 5)
```

### Integrate datasets
```{r Integrate datasets}
## Run Liger iNMF integration

## Normalize data, select variable genes, and scale data
plan('sequential')
so <- rliger::normalize(so) 
so <- selectGenes(so, datasetVar='sample', thresh=0.1)
so <- scaleNotCenter(so, datasetVar='sample')

## Run integrated iNMF and align factors
gc()
so <- runINMF(so, datasetVar='sample', k = 20, nCores=20, lambda=5)
so <- alignFactors(so, method = "centroidAlign")
```

## UMAP dimensional reduction and broad clustering
```{r UMAP dimensional reduction and broad clustering}
## Dimensional reduction and simple clustering (B cells and T cells)

## UMAP dimensional reduction
plan('multicore')
so <- RunUMAP(so, reduction = "inmfNorm", dims = 1:20, min.dist=0.01,
              seed.use=50)


## Create neighborhood graph from normalized inmf embedding
so <- FindNeighbors(so, reduction='inmfNorm', dims=1:20, l2.norm=T, 
                    graph.name=c('RNA_nn','RNA_snn'))

## Determine optimal clustering
# Iteratively assess different clustering resolutions 
plan('sequential')
oc <- DetermineOptimalClusters(so, graph.name='RNA_snn', 
                               resolution=seq(0.05,1.0,0.05))
# Plot WSS and silhouette-widths for each resolution
oc$plot


## Cluster cells 
so <- FindClusters(so, graph.name='RNA_snn', algorithm=4, resolution=0.4, 
                   cluster.name='clust', random.seed=42)
DimPlot(so, label=T, group.by='clust')


## Save data
saveRDS(so, 'so_proc.rds')
```

## Visualize and filter out B cells
```{r Visualize and filter out B cells}
## Visualize and filter out B cells

## Add umap embedding coordinates to the metadata of the seurat object
so$x <- so$umap@cell.embeddings[,1]
so$y <- so$umap@cell.embeddings[,2]


## Plot UMAP of all cells (For schematic)
p <- ggplot(so@meta.data, aes(x=x,y=y)) +
    geom_point_rast(size=0.1, color='gray50') +
    geom_point_rast(data= so@meta.data[which(so$ct !='Bcell'),], 
                    size=0.01, color='steelblue') +
    theme_blank() + theme(legend.position='none')
ggsave(plot=p, 'Figures/umap_allcells.pdf', height=2.5*1, width=2.5*1,
       scale=1.5)


## Visualize bulk B and T cell populations on UMAP
l <- data.frame(x=c(-12.5,-5,-5), y=c(-1,-1,8))
l <- geom_line(data=l, aes(x=x,y=y), linetype='dashed')
p <- list(FeaturePlot(so,'MS4A1', raster=T) + l, 
          FeaturePlot(so,'CD3E', raster=T) + l)
p <- cowplot::plot_grid(plotlist=p, align='v', axis='lbt', nrow=2)
ggsave(plot=p, 'Figures/_BcellRemoval_step1.pdf', height=2.5*2, width=2.9*1,
       scale=1.5)


## Exclude bulk B cell population
cells2keep <- colnames(so)[which(!(so$x < -5 & so$y > -1))]
so <- subset(so, cells=cells2keep)
so <- subset(so, subset=(ct != 'Bcell'))
```

## Cluster cell subsets
```{r Cluster cell subsets}

## Normalize and scale data
plan('multicore')
so <- NormalizeData(so)
plan('sequential')
so <- ScaleData(so)

## load library
library(monocle3)

## convert so to cds to use Monocle 3 clustering
cds <- SeuratWrappers::as.cell_data_set(so, assay = 'RNA', 
                                        reductions = Reductions(so))
rowData(cds)$gene_short_name <- row.names(rowData(cds))


## Preprocess data with Monocle 3 (none of this is used, but Monocle 3 reqs it)
cds <- preprocess_cds(cds, method = 'PCA', norm_method = 'none', 
                      nn_control=list(nn.cores=20))


## Cluster and visualize with Monocle 3 
cds <- cluster_cells(cds, reduction_method='UMAP')
plot_cells(cds, show_trajectory_graph=F) 


## Add cluster numbers to SO
so$c2=clusters(cds)
Idents(so) <- so$c2
```

## Classify cell populations
```{r Classify cell populations}

## Create list of pathway gene ontologies from MSigDB (Reactome and Hallmark)
gmt <- list(fgsea::gmtPathways(
    'SCT/00_ReferenceInfo/c2.reactome_pathways.v7.4.symbols.gmt.txt'),
            fgsea::gmtPathways(
                'SCT/00_ReferenceInfo/MSigDB_Hallmark_2020.txt')) %>% 
    unlist(.,F)
# Filter pathway genes by those with count data
gmt <- lapply(gmt, function(x) intersect(x, rownames(so)))


## Create gene lists with manual cell-type ontologies and pathway ontologies
glist2 <- list(Effector=c('TBX21','KLRG1','ZEB2','GZMA','GZMB','FASLG'),
               TRM=c('ZNF683','CXCR3','CXCR6'),
               Exh=c('TIGIT','HAVCR2','TOX','ENTPD1'),
               MEM=c('KLF2','SELL','EOMES','S1PR1','IL15RA'), #,
               TCM_CD4=c('CCR7','TCF7','CD4'),
               TPEX=c('SLAMF6','TOX'),
               BCell=c('CD19','MS4A1'),
               TCell=c('CD3D','CD3E'),
               TH1=c('CD4','ZEB2','IFNG','TBX21','TNF','GATA3','IL2RB'),
               TC=c('CD8A','CD8B'),
               TH=c('CD4'),
               Cytotoxic=c('GZMA','GZMB','PRF1','FASLG'),
               Cytokines=c('TNF','IFNG'),
               Prolif=c('TOP2A','MKI67'),
               OXPHOS=gmt[['Oxidative Phosphorylation']],
               Glycol=gmt[['Glycolysis']]
)


## Create gene module scores of glist2 celltype and metabolic gene sets
# Create a matrix of scores, update colnames, then add to SO metadata
plan('sequential')
DefaultAssay(so) <- 'RNA'
temp <- AddModuleScore(so, features=glist2, name='glist2_', seed=1)@meta.data
temp <- temp[,grep('glist2',colnames(temp))]
colnames(temp) <- names(glist2)
for(x in names(glist2)) so@meta.data[,x] <- temp[,x]
rm(temp)

## Rename clusters
so$ct <- c('TE_1','TE_2','TEX','TE_3','TRM_1',
           'TMEM','TRM_2','TPEX','TH1','TE_4',
           'TE_5','Act.T','Bcell','TH1_TCM')[as.integer(so$c2)]


## Create grouped subsets
so$groups <- str_split_i(so$ct,'_',1)
so$groups[which(so$groups =='TH1')] <- 'CD4T'
so$groups[which(so$groups =='TEX')] <- 'TPEX/TEX'
so$groups[which(so$groups =='TPEX')] <- 'TPEX/TEX'
Idents(so) <- so$ct
Idents(so) <- factor(so$ct)
```

### Exclude ambiguous cells
```{r Exclude ambiguous cells}
## Scattered B cells were found within Tcell-majority clusters. These were 
##    termed "ambiguous" and removed below.

## label and exclude ambiguous
so$ambiguous <- (colSums(so$RNA$counts[c(grep('^IGH', rownames(so), value=T)
                                         ,'MS4A1', 'CD19', 'CD22'),]) > 0.01)
# Make sure the B cell cluster is not labeled as ambiguous
so$ambiguous[which(so$ct=='Bcell')] <- F


## Exclude ambigous cells
so <- subset(so, subset=(ambiguous == F ))

```

### Visualize marker genes and signatures across cell subsets
```{r Visualize marker genes and signatures across cell subsets}

## Create a filtered data set of just cytotoxic T cells (so2)
celllist1 <- c('Bcell','TH1','CD4_TCM', 'TMEM','TRM_1','TRM_2',
  paste0('TE_',1:5),'Act.T','TEX','TPEX')
celllist2 <- c('TMEM','TRM_1','TRM_2',paste0('TE_',1:5),
               'Act.T','TPEX','TEX')
so2 <- subset(so, subset=(ct %in% celllist2))


## Make custom dot plots
p <- list(
    # Plot of major cell populations
    customDotPlot(so, c('BCell','TCell','TC','TH','TCM_CD4'),
                  group.order=rev(celllist1), scale.max=50),
    # Plot of CD8 T cell subsets and related signatures
    customDotPlot(so2, 
                  c('MEM','TRM','Effector','TPEX','Exh','Prolif'),
                  rev(cells2display)),
    # Plot of CD8 T cell subsets and phenotype signatures
    customDotPlot(so2, 
                  c('Glycol','OXPHOS','Cytotoxic','Cytokines'),
                  rev(cells2display)))
p <- cowplot::plot_grid(plotlist=p, align='h', axis='lbt', nrow=1)
ggsave(plot=p, 'Figures/dotplot_celltypes_v0.4.pdf', 
       height=3*1, width=1.5*3,scale=1.5)


## Make a dotplot of the genes used for the subset classification signatures
# Create marker lists
glist1 <- list(BCell=c('CD19','MS4A1'),
              TCell=c('CD3D','CD3E'),
              TC=c('CD8A','CD8B'),
              TH=c('CD4'),
              TCM_CD4=c('CCR7','TCF7')
) %>% unlist() %>% unique()
glist2 <- list(MEM=c('KLF2','SELL','EOMES','S1PR1','IL15RA'),
              TRM=c('ZNF683','CXCR3','CXCR6'),
              Effector=c('TBX21','KLRG1','ZEB2','GZMA','GZMB','FASLG'),
              TPEX=c('SLAMF6','TOX'),
              Exh=c('TIGIT','HAVCR2','ENTPD1'),
              Prolif=c('TOP2A','MKI67'),
              Cytokines=c('TNF','IFNG'),
              Glycol=c('GAPDH','PGAM1','LDHA'),
              OXPHOS=c('MDH1','MDH2','IDH2')
) %>% unlist() %>% unique()
# create plot lists
plist <- list(
    # major cell subset dot plot
    customDotPlot(so, glist1,rev(celllist1)),
    # CD8-specific dot plot
    customDotPlot(so2, glist2,rev(celllist2))
) 
# Export plots
p <- cowplot::plot_grid(plotlist=plist, align='h', axis='l', nrow=1, rel_widths=c(0.4,1))
ggsave(plot=p, 'Figures/dotplot_celltypes_genes_v0.5.pdf', 
       height=2.8, width=9,scale=1.5)

## remove the filtered SO 
rm(so2)
```

### UMAP plots of processed data
```{r UMAP plots of processed data}

## Remove B cells
so <- subset(so, subset=(ct != 'Bcell'))

## Save data
saveRDS(so, 'Datasets/so_proc1.rds')

## Plot UMAPs
# find center coordinates of each cluster
centroids <- reframe(so@meta.data, .by=c('ct'), 
                     label=ct[1], x=mean(x), y=mean(y))
# Make plot of complete data
p1 <- ggplot(so@meta.data, aes(x=x,y=y)) +
    ggrastr::geom_point_rast(size=2,color='gray25') +
    ggrastr::geom_point_rast(size=0.75, aes(color=ct)) + theme_void() +
    ggrepel::geom_text_repel(data=centroids, aes(x=x,y=y,label=label),
                             min.segment.length = 0) +
    theme(legend.position='none')
# Make separate plots for each genotype
cols <- structure(ggcolor(14), names=sort(unique(so$ct)))
p2 <- lapply(c('GTC','NR4A1','SOCS3','DKO'), function(celltype){
    d <- so@meta.data[which(so$sample==celltype),]
    centroids <- reframe(d, .by=c('ct'), label=ct[1], x=mean(x), y=mean(y))
    ggplot(d, aes(x=x,y=y)) + theme_void() +
        ggrastr::geom_point_rast(data=so@meta.data,size=1,color='gray25') +
        ggrastr::geom_point_rast(data=so@meta.data,size=0.25, color='gray98') + 
        ggrastr::geom_point_rast(size=1,color='gray25') +
        ggrastr::geom_point_rast(size=0.25, aes(color=ct)) + 
        labs(title=celltype) +
        scale_color_manual(values=cols) + 
        theme(legend.position='none', title=element_text(vjust=0.25))
})
p2 <- cowplot::plot_grid(plotlist=p2, align='hv', axis='lb')
# Combine the plots and export
p <- cowplot::plot_grid(p1,p2, align='hv', axis='lb')
ggsave(plot=p, 'Figures/umap_celltypes_v0.1.pdf', height=3*1, width=3*2)

```

### Compare distributions of cell subsets
```{r Compare distributions of cell subsets}

## Make stacked bar plots to compare celltype %s across genotypes
# Make dataframe 'd' with % of cells belonging to each celltype, separately 
#   for each genotype
d <- so@meta.data[,c('ct', 'groups', 'sample')]
d <- mutate(d, .by=c('sample'), total=length(sample))
d <- reframe(d, .by=c('ct','sample'), pct=length(ct)/total[1])
# Set the ploting order for each celltype and genotype
d$sample <- factor(d$sample, levels=rev(c('GTC','NR4A1','SOCS3','DKO')))
celltypes <- c("TMEM","TH1","CD4_TCM","TRM_1","TRM_2",
               "TE_1","TE_2","TE_3","TE_4","TE_5",
               "Act.T","TPEX","TEX")
d$ct <- factor(d$ct, levels=rev(celltypes))
# Make stacked bar plot 
cols <- structure(ggcolor(14), names=sort(unique(so$ct)))
p1 <- ggplot(d, aes(y=sample, x=pct, fill=ct)) + 
    geom_bar(stat='identity', position='stack', color='gray20') + theme_classic() +
    labs(y=NULL, x=('Proportion of T cells')) + 
    theme(legend.key.spacing=unit(0.5,'lines'), 
          legend.key.size = unit(0.9, 'lines'),
          legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
          legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
          legend.title=element_text(size=10))

## Make non-stacked bar plots to compare celltype %s across genotypes
# Make dataframe 'd' with % of cells belonging to each celltype, separately 
#   for each genotype
# Note: removed 'CD4_TCM', because it is only represented by GTC genotype and 
#   thus, cannot be compared to other genotypes.
d <- so@meta.data[which(so$ct != 'CD4_TCM'),c('ct', 'groups', 'sample')]
d <- mutate(d, .by=c('sample'), total=length(sample))
d <- reframe(d, .by=c('ct','sample'), pct=length(groups)/total[1])
d$ct <- factor(d$ct, levels=celltypes)
# Set the ploting order for each celltype and genotype
d$sample <- factor(d$sample, levels=(c('GTC','NR4A1','SOCS3','DKO')))
cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                  names=c('GTC','SOCS3','NR4A1','DKO'))
# Make bar plot
p2 <- ggplot(d, aes(x=ct, y=pct, fill=sample)) + 
    geom_bar(stat='identity', position='dodge', color='gray20') + theme_classic() +
    scale_fill_manual(values=cols) +
    labs(x=NULL, y=str_wrap('Proportion of T cells')) + #coord_flip()+
    theme(legend.key.spacing=unit(0.5,'lines'), #legend.position='bottom', 
          legend.key.size = unit(0.9, 'lines'),
          legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
          legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
          #axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_text(size=10))
# Combine plots and export
p <- cowplot::plot_grid(p1,p2, ncol=1, align='hv', axis='lbt', rel_heights=c(0.4,0.6))
ggsave(plot=p, 'Figures/barplot_subset_pct_v0.4.pdf', height=1.6*2, width=8)


## Make stacked bar plots to compare celltype group %s across genotypes
# Make dataframe 'd' with % of cells belonging to each group, separately 
#   for each genotype
d <- so@meta.data[ ,c('ct', 'groups', 'sample')]
d <- mutate(d, .by=c('sample'), total=length(sample))
d <- reframe(d, .by=c('groups','sample'), pct=length(groups)/total[1])
# Set the ploting order for each celltype and genotype
d$sample <- factor(d$sample, levels=rev(c('GTC','NR4A1','SOCS3','DKO')))
d$groups <- factor(d$groups, levels=rev(c('TMEM','TRM','TE','Act.T','TPEX/TEX','TH1')))
# Make stacked bar plot
p1 <- ggplot(d, aes(y=sample, x=pct, fill=groups)) + 
    geom_bar(stat='identity', position='stack', color='gray20') + theme_classic() +
    #scale_fill_manual(values=cols) +
    labs(y=NULL, x=('Proportion of T cells')) + #coord_flip()+
    #guides(fill=guide_legend(ncol=2)) +
    theme(legend.key.spacing=unit(0.5,'lines'), #legend.position='bottom', 
          legend.key.size = unit(0.9, 'lines'),
          legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
          legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
          #axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_text(size=10))


## Make non-stacked bar plots to compare celltype %s across genotypes
# Make dataframe 'd' with % of cells belonging to each celltype, separately 
#   for each genotype
d <- so@meta.data[,c('ct', 'groups', 'sample')]
d <- mutate(d, .by=c('sample'), total=length(sample))
d <- reframe(d, .by=c('groups','sample'), pct=length(groups)/total[1])
# Set the ploting order for each celltype and genotype
d$groups <- factor(d$groups, levels=c('TMEM','TRM','TE','Act.T','TPEX/TEX','TH1'))
d$sample <- factor(d$sample, levels=(c('GTC','NR4A1','SOCS3','DKO')))
cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                  names=c('GTC','SOCS3','NR4A1','DKO'))
# Make bar plot
p2 <- ggplot(d, aes(x=groups, y=pct, fill=sample)) + 
    geom_bar(stat='identity', position='dodge', color='gray20') + theme_classic() +
    scale_fill_manual(values=cols) +
    labs(x=NULL, y=str_wrap('Proportion of T cells')) + #coord_flip()+
    theme(legend.key.spacing=unit(0.5,'lines'), #legend.position='bottom', 
          legend.key.size = unit(0.9, 'lines'),
          legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
          legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
          #axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_text(size=10))
# Combine plots and export
p <- cowplot::plot_grid(p1,p2, ncol=1, align='hv', axis='lbt', rel_heights=c(0.4,0.6))
ggsave(plot=p, 'Figures/barplot_subset_pct_groups_v0.2c.pdf', height=2*2, width=8)



## Generate stats for celltype proportions 
library(rstatix)
# generate and arrange data: # cells in and out of each subset for each genotype
d <- so@meta.data[ ,c( 'ct', 'sample')]
d <- reframe(d, .by=c('ct','sample'), x=length(ct))
d <- mutate(d, .by='sample', total=sum(x))
d$y <- d$total - d$x
rownames(d) <- paste0(d$sample,'_',d$ct)
# Run pairwise fisher exact tests on counts
r <- do.call(rbind, lapply(celltypes, function(ct){
    if(ct=='CD4_TCM') return(NULL)
    r <- pairwise_fisher_test(d[which(d$ct==ct), c('x','y')])
    return(data.frame(ct=ct, r[,1:4]))
}))
# Adjust p values for multiple testing
r$p.adj <- p.adjust(r$p)
# Add column to sort results in order of the bar plot
comps <- c('GTC-DKO','NR4A1-DKO','SOCS3-DKO','GTC-NR4A1','GTC-SOCS3','NR4A1-SOCS3')
geno <- c('GTC','NR4A1','SOCS3','DKO')
r$comp <- apply(r[,2:3],1, function(x) {
    x <- str_split_i(x,'_',1)
    paste0(as.character(sort(factor(x,levels=geno))),collapse='-')
})
r$comp <- factor(r$comp, levels=comps)
r <- r[order(r$ct, r$comp),]
# Add column to show astrices for signif.
r$signif <- case_when(r$p.adj < 1e-4 ~ '****',
                      r$p.adj < 1e-3 ~ '***',
                      r$p.adj < 1e-2 ~ '**',
                      r$p.adj < 0.05 ~ '*',
                      r$p.adj >= 0.05 ~'ns')
# Export table
write.table(r, 'Data/subset_pct_comp.txt', sep='\t', quote=F, row.names=F)


## Generate stats for Grouped celltype proportions 
library(rstatix)
# generate and arrange data: # cells in and out of each subset for each genotype
d <- so@meta.data[ ,c( 'groups', 'sample')]
d <- reframe(d, .by=c('groups','sample'), x=length(groups))
d <- mutate(d, .by='sample', total=sum(x))
d$y <- d$total - d$x
d$p <- d$x/d$total
rownames(d) <- paste0(d$sample,'_',d$groups)
# Run pairwise fisher exact tests on counts
r <- do.call(rbind, lapply(unique(d$groups), function(gp){
    r <- pairwise_fisher_test(d[which(d$groups==gp), c('x','y')])
    return(data.frame(ct=gp, r[,1:4]))
}))
# Adjust p values for multiple testing
r$p.adj <- p.adjust(r$p)
# Add column to sort results in order of the bar plot
comps <- c('GTC-DKO','NR4A1-DKO','SOCS3-DKO','GTC-NR4A1','GTC-SOCS3','NR4A1-SOCS3')
geno <- c('GTC','NR4A1','SOCS3','DKO')
r$comp <- apply(r[,2:3],1, function(x) {
    x <- str_split_i(x,'_',1)
    paste0(as.character(sort(factor(x,levels=geno))),collapse='-')
})
r$comp <- factor(r$comp, levels=comps)
r <- r[order(r$ct, r$comp),]
# Add column to show astrices for signif.
r$signif <- case_when(r$p.adj < 1e-4 ~ '****',
                      r$p.adj < 1e-3 ~ '***',
                      r$p.adj < 1e-2 ~ '**',
                      r$p.adj < 0.05 ~ '*',
                      r$p.adj >= 0.05 ~'ns')
# Export table
write.table(r, 'Data/subsetGroups_pct_comp.txt', sep='\t', quote=F, row.names=F)
```
