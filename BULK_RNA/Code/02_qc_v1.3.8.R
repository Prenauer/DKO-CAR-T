library(edgeR)
library(stringr)
library(dplyr)
require(ggrepel)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(DESeq2)
library(cowplot)
options(stringsAsFactors = F)
setwd('dkot_rna')


## Accessory function for plotting
ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]


###################

## Setup data for DE
data <- catchKallisto(list.dirs('Data/kallisto')[-1])
colnames(data$counts) <- apply(stringr::str_split(basename(colnames(data$counts)),'_',simplify=T)[,1:3], 1, function(x) paste0(x, collapse='_'))
meta <- data.frame(colnames(data$counts), stringr::str_split(colnames(data$counts), '_', simplify=T), row.names=colnames(data$counts))
colnames(meta) <- c('sample','genotype','time','replicate')
meta$genotype <- c('DKO','GTC','NR4A1','SOCS3')[factor(meta$genotype)]
meta$genotype <- factor(meta$genotype, levels=c('GTC','NR4A1','SOCS3','DKO'))
meta$time <- factor(meta$time, levels=c('0','6','24','72'))
meta <- meta[with(meta, order(time, genotype, replicate)),]
write.table(data.frame(ENST=rownames(data$counts), data$counts[,meta$sample]), 'Data/02a_dkotrna_txcounts.txt', sep = '\t', col.names = T, row.names = F, quote = F)
write.table(meta, 'Data/02b_dkotrna_metadata_table.txt', sep = '\t', col.names = T, row.names = T, quote = F)


## setup gene annotations
## annotate genes
gtf <- rtracklayer::import('Downloads/Homo_sapiens.GRCh38.96.gtf')
gtf <- data.frame(gtf)
gtf <- gtf[!is.na(gtf$transcript_id),]
gtf$transcript_id <- paste0(gtf$transcript_id,'.', gtf$transcript_version)
gtf$gene_id <- paste0(gtf$gene_id,'.', gtf$gene_version)
gtf <- gtf[,c('gene_id','gene_name','transcript_id','transcript_name','transcript_support_level')] %>% unique()
# add missing tx to the file
gtf <- rbind(gtf, data.frame(gene_id=setdiff(data$ENST, gtf$transcript_id), 
                             gene_name=NA, 
                             transcript_id=setdiff(data$ENST, gtf$transcript_id),
                             transcript_name=NA,
                             transcript_support_level=NA,
                             row.names=setdiff(data$ENST, gtf$transcript_id)))
gtf$transcript_support_level <- stringr::str_split(gtf$transcript_support_level,' ', simplify=T)[,1]
rownames(gtf) <- gtf$transcript_id
write.table(gtf[data$ENST,], 'Reference/gencode_GRCh38.96_tx_info.txt', sep = '\t', col.names = T, row.names = T, quote = F)


###########################

## load prepped data
data <- read.delim('02a_dkotrna_txcounts.txt', row.names=1)
meta <- read.delim('02b_dkotrna_metadata_table.txt')
rownames(meta) <- meta$sample
meta$sample <- factor(meta$sample, levels=meta$sample)
meta$genotype <- factor(meta$genotype, levels=c('GTC','NR4A1','SOCS3','DKO'))
meta$time <- factor(meta$time, levels=c('0','6','24','72'))
gtf <- read.delim('Reference/gencode_GRCh38.96_tx_info.txt')
tx2use <- gtf$transcript_id[which(gtf$transcript_support_level < 4)]


## Correlation plot of all samples
lcpm <- DESeq2::varianceStabilizingTransformation(round(as.matrix(data[filterByExpr(data),])), fitType='local')
lcpm <- merge(gtf[,3:4], lcpm, by.x='transcript_id', by.y=0)
cor.data <- pcaPP::cor.fk(lcpm[, 3:ncol(lcpm)])
hm.filename <-  paste0('Figures/02a_hm_allgenes_allConditions.pdf')
hm <- pheatmap::pheatmap(cor.data, scale='none', cluster_cols=T, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
                         annotation_row = meta[,2:3], annotation_col = meta[,2:3], 
                         show_rownames=T, treeheight_col=15, treeheight_row=15, clustering_method='average', main='Sample correlation (Kendall)',
                         angle_col='90', cutree_rows=5, cutree_cols=5, filename=hm.filename, height=7, width=8)


## Correlation plot for each timepoint
pdf('Figures/02b_hm_corrplot_allgenes_byTime.pdf', height=3.2, width=4.2)
for(timepoint in c(0,6,24,72)) {
    d <- as.matrix(data[filterByExpr(data), as.character(meta$sample[which(meta$time==timepoint)])])
    d <- DESeq2::varianceStabilizingTransformation(round(d), fitType='local')
    cor.data <- pcaPP::cor.fk(d)
    pheatmap::pheatmap(cor.data, scale='none', cluster_cols=T, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
                       show_rownames=T, treeheight_col=15, treeheight_row=15, clustering_method='average', main=paste0(timepoint, ' hrs activation'),
                       angle_col='90', cutree_rows=1, cutree_cols=1, annotation_row = meta[,2, drop=F], annotation_col = meta[,2, drop=F])
}
dev.off()



## MDS Plot of all samples
pca <- prcomp(scale(t(lcpm[, 3:ncol(lcpm)])))$x[,1:2] 
pca <- data.frame(pca, meta[, 1:3])
p <- ggplot(pca, aes(x=PC1, y=PC2, color=genotype)) + geom_point() + 
    ggrepel::geom_text_repel(aes(label=sample),size=2.25, min.segment.length=0, nudge_y=0.5, segment.size=0.2,segment.alpha=0.4,max.overlaps=50) +
    theme_test() + theme(legend.position='none') + labs(title='PCA of all samples') 
ggsave(plot=p, 'Figures/02c_mds_allSamples.pdf', height=3, width=3.25)
p <- ggplot(pca, aes(x=PC1, y=PC2, color=genotype)) + geom_point() + 
    ggforce::geom_mark_hull(aes(group=time), color='gray50', concavity=2, expand = unit(2.5, "mm")) + 
    scale_color_manual(values=structure(ggcolor(4), names=c('DKO','SOCS3','NR4A1','GTC'))) +
    theme_test() + theme(legend.position='none') + labs(title='PCA of all samples') 
ggsave(plot=p, 'Figures/02c_mds_allSamples_hull.pdf', height=3.5, width=4)



## MDS Plot for each timepoint
p <- lapply(c(0,6,24,72), function(timepoint){
    pca <- prcomp(scale(t(lcpm[, as.character(meta$sample[which(meta$time==timepoint)])])))$x[,1:2] 
    pca <- data.frame(pca, meta[which(meta$time==timepoint), 1:3])
    ggplot(pca, aes(x=PC1, y=PC2, color=genotype)) + geom_point() + 
        ggrepel::geom_text_repel(aes(label=sample),size=3, min.segment.length=0, nudge_y=0.5,segment.alpha=0.4, segment.size=0.2) +
        theme_test() + theme(legend.position='none') + labs(title=paste0(timepoint, ' hrs activation')) %>% return()
})
p <- plot_grid(plotlist=p, align='h', ncol=2)
ggsave(plot=p, 'Figures/02d_mds_samples_byTime.pdf', height=4, width=4.5)
p <- lapply(c(0,6,24,72), function(timepoint){
    pca <- prcomp(scale(t(lcpm[, as.character(meta$sample[which(meta$time==timepoint)])])))$x[,1:2] 
    pca <- data.frame(pca, meta[which(meta$time==timepoint), 1:3])
    ggplot(pca, aes(x=PC1, y=PC2, color=genotype)) + geom_point() + 
        ggforce::geom_mark_hull(aes(group=genotype), expand = unit(2.5, "mm")) + 
        scale_color_manual(values=structure(ggcolor(4), names=c('DKO','SOCS3','NR4A1','GTC'))) +
        theme_test() + theme(legend.position='none') + labs(title=paste0(timepoint, ' hrs activation')) %>% return()
})
p <- plot_grid(plotlist=p, align='h', ncol=2)
ggsave(plot=p, 'Figures/02d_mds_samples_byTime_hull.pdf', height=3.5, width=4.5)


