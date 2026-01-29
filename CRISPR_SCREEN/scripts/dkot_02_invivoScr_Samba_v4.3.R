library(edgeR)
library(stringr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(reshape2)
library(factoextra)
library(cluster)
options(ggrastr.default.dpi=750)
source('Code/samba_v2.0.R')


ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
addSmallLegend <- function(size = NULL, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    if(!is.null(size)){
        pointSize = pointSize * size
        textSize = textSize * size
        spaceLegend = spaceLegend * size
    }
    list(guides(shape = guide_legend(override.aes = list(size = pointSize)),
                color = guide_legend(override.aes = list(size = pointSize))),
         theme(legend.title = element_text(size = textSize), 
               legend.text  = element_text(size = textSize),
               legend.key.size = unit(spaceLegend, "lines")))
}


## Load data
rc <- read.delim('Data/counts/dkot_invivoScrTumor_v2m2.txt')
rc[is.na(rc)] <- 0


## setup design matrix
tumor <- grepl('T',colnames(rc))[3:ncol(rc)] %>% as.integer()
cell <- as.integer(!grepl('Plasmid',colnames(rc))[3:ncol(rc)]) 
donor1 <- as.integer(substr(colnames(rc)[3:ncol(rc)],2,2) == '3')
donor2 <- as.integer(substr(colnames(rc)[3:ncol(rc)],2,2) == '5')
donor3 <- as.integer(substr(colnames(rc)[3:ncol(rc)],2,2) == '6')
design <- model.matrix(~ cell + donor1 + donor2 + tumor)



## Complete analysis
dge <- Preprocess_Samba(data = rc, design = design, group = tumor, min.guides = 1, normalization.method = 'uq2')
sgRes <- Analyze_Samba_Guides(dge = dge, coefficient ='tumor', method ='QLF', file.prefix = 'Data/dko_invivoScr')
sgRes <- sgRes[with(sgRes, order(PValue, decreasing = F)),]
sgRes[which(sgRes$Gene == 'GTC'),'Gene'] <- 'NTC'
geneRes <- Analyze_Samba_Genes(sgRes = sgRes, ntc.as.null.dist = F, score.method ='GeneScore', file.prefix = 'Data/dko_invivoScr')
saveRDS(dge, file = paste0('Data/dko_invivoScr_dge.rds'))


## Save normalized data
df <- merge(rc[,1:2], dge$counts, by.x = 1, by.y = 0, all.x = F, all.y = T)
write.table(df, 'Data/dkot_invivo_counts_uq2-norm.txt', sep = '\t', row.names = F, quote = F)
df <- cbind(rc[,1:2], log2(1+apply(rc[,3:ncol(rc)],2, function(x) 1e6*x/sum(x))))
write.table(df, 'Data/dkot_invivo_counts_log2cpm-norm.txt', sep = '\t', row.names = F, quote = F)


## Random plot of SAMBA hits
geneRes <- read.delim('Data/dko_invivoScr_GeneLevelResults.txt')
geneRes <- geneRes[order(geneRes$pval_pos, decreasing=F),]
geneRes <- geneRes[which(geneRes$Gene != 'GTC'),]
set.seed(1)
d <- data.frame(gene = geneRes$Gene,
                logQ = -log10(geneRes$qval_pos),
                z = scale(geneRes$score_pos),
                nguides = geneRes$n_fdr_guides_pos,
                rank = order(geneRes$score_pos, decreasing = T),
                geneorder = sample(1:nrow(geneRes), nrow(geneRes), T))
d$gene <- factor(d$gene)
top_sko <- d[grep('_',d$gene, invert=T),]
top_sko <- top_sko[order(top_sko$rank),] %>% head(2)
z.thresh <- 2.84 #min(d$z[which(d$logQ > -log10(0.05))] )
d.label <- rbind(d[which(abs(d$z) > z.thresh),], top_sko) %>% unique()
d.label$gene <- factor(d.label$gene, levels = unique(d.label$gene))
pt.scaling.factor <- 0.1 
p <- ggplot(d, aes(x = geneorder, y = z), alpha = 0.25) + 
    rasterize(ggpointdensity::geom_pointdensity(show.legend = F, adjust=10, aes(size = nguides*pt.scaling.factor))) + 
    scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
    scale_size_continuous(range = c(0.2,2.5), name='# FDR guides') + 
    geom_hline(yintercept = c(z.thresh,-z.thresh), linetype = 'dashed', color = 'gray30', alpha=0.4, linewidth = 0.4) +

    ggnewscale::new_scale_color() + 
    geom_point_rast(data = d.label, aes(size = 1.5*nguides*pt.scaling.factor), color = 'gray10') + 
    geom_point_rast(data = d.label, aes(size = nguides*pt.scaling.factor), color = 'white') + 
    geom_point_rast(data = d.label[which(d.label$z >  z.thresh),], aes(size = nguides*pt.scaling.factor), color = 'steelblue', alpha=1)+ 
    geom_point_rast(data = d.label[which(d.label$z < -z.thresh),], aes(size = nguides*pt.scaling.factor), color = 'firebrick', alpha=1) + 
    ggrepel::geom_text_repel(data = head(d.label,10), segment.alpha = 0.5, aes(label = gene), #fontface='italic', 
                             size = 2.8, max.iter = 1000000, nudge_y=0.5, nudge_x=1,
                             min.segment.length = 0, segment.size = 0.1, max.overlaps = 10) +
    ggrepel::geom_text_repel(data = tail(d.label,6), segment.alpha = 0.5, aes(label = gene), #fontface='italic', 
                             size = 2.8, max.iter = 1000000, nudge_y=-1, nudge_x=1,
                             min.segment.length = 0, segment.size = 0.1, max.overlaps = 10) +
    ylim(c(-5,6)) +
    labs(x = 'Gene number', y = 'Gene z-score (Samba)') + 
    cowplot::theme_cowplot() + addSmallLegend(3)
ggsave(plot = p, filename = 'Figs/dkot_invivo_randplot_v4.3.2.pdf', width = 5, height = 3.2)


## MDS plots
dge <- readRDS(file = paste0('Data/dkot_invivoScr_dge.rds'))
lcpm <- edgeR::cpm(dge, log = T) %>% scale() %>% data.frame() #%>% scale() 
lcpm.batchremoved <- limma::removeBatchEffect(lcpm, batch = tumor)
mds <- limma::plotMDS(lcpm, top = 10000,gene.selection = 'common')
axis.labels <- paste0(c('Dim 1 (', 'Dim 2 ('),round(mds$var.explained[1:2] * 100, 1), c('% variance)', '% variance)'))
mds <- data.frame(Dim1 = mds$x, Dim2 = mds$y, row.names = colnames(lcpm))
set.seed(42)
km.res <- kmeans(mds, 3, nstart = 25)
p1 <- fviz_cluster(km.res, data = mds, shape = 19, xlab = axis.labels[1], ylab = axis.labels[2],
                   palette = c("#00AFBB", "#E7B800", "#FC4E07",'gray50'),
                   ggtheme = theme_test(), show.clust.cent=F, repel=T, labelsize=8, ylim=c(-1.5,2), xlim=c(-1.5,2.25),
                   main = 'Unadjusted MDS')
data.output <- data.frame(Plot='Unadjusted MDS', Sample=rownames(mds), mds, km.res$cluster)
mds <- limma::plotMDS(lcpm.batchremoved, top =10000, gene.selection = 'common')
axis.labels <- paste0(c('Dim 1 (', 'Dim 2 ('),round(mds$var.explained[1:2] * 100, 1), c('% variance)', '% variance)'))
mds <- data.frame(Dim1 = mds$x, Dim2 = mds$y, row.names = colnames(lcpm))
set.seed(42)
km.res <- kmeans(mds, 3, nstart = 25)
p2 <- fviz_cluster(km.res, data = mds, shape = 19, xlab = axis.labels[1], ylab = axis.labels[2], ellipse.alpha=0.1,
                   palette = c("#00AFBB", "#E7B800", "#FC4E07"), #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                   ggtheme = theme_test(), show.clust.cent=F, repel=T, labelsize=8, ylim=c(-3,1.5), xlim=c(-1.25,2.25),
                   main = 'Coef-adjusted MDS')
p1 | p2
ggsave(plot = p1 | p2, filename = 'Figs/qc_invivo_mds_v2.pdf', width = 7, height = 3)
data.output <- rbind(data.output, data.frame(Plot='Adjusted MDS', Sample=rownames(mds), mds, km.res$cluster))
write.table(data.output, 'Data/dko_invivoScr_mds_data.txt', sep = '\t', quote = F, row.names = F)


##################

## Plot top sgRNAs by sample

# Get top gene-pairs
geneRes <- read.delim('Data/dko_invivoScr_GeneLevelResults.txt')
geneRes <- geneRes[order(geneRes$score_pos, decreasing=T),]
tophits <- geneRes$Gene[order(geneRes$pval_pos)] %>% head(6)


# Get top sgRNAs of gene-pairs
sgRes <- read.delim('Data/dko_invivoScr_GuideLevelResults.txt')
sgRes <- sgRes[order(sgRes$logFC, decreasing=T),]
top4guide <- slice_head(sgRes[which(sgRes$Gene %in% tophits),], by='Gene', n=4)[,1]


# Get data
dge <- readRDS(file = paste0('Data/dkot_invivoScr_dge.rds'))
rc <- read.delim('Data/dkot_invivoScrTumor_v2m2.txt')
rc <- rc[which(rc$sgRNA %in% sgRes$sgRNA),]

# z-score data
lcpm <- edgeR::cpm(dge, log = T) 
lcpm <- merge(rc[,1:2], lcpm, by.x = 1, by.y = 0, all.x = F, all.y = T)
gtc.medians <- colMeds(lcpm[which(lcpm$Gene == 'GTC'), 3:ncol(lcpm)])
z <- cbind(lcpm[,1:2,], data.frame(scale(lcpm[, 3:ncol(lcpm)], center=gtc.medians, scale=T)))
gtc <- data.frame(sgRNA='GTC', Gene='GTC', matrix(colMeds(z[which(z$Gene == 'GTC'), 3:ncol(z)]), nrow=1))
colnames(gtc) <- colnames(z)
# reshape data
d <- rbind(z[which(z$sgRNA %in% top4guide), ], gtc)
d <- reshape2::melt(d)
colnames(d) <- c('sgRNA','Gene','Sample','Value')
# order the sgRNA, Genes, and samples
d$Gene <- factor(d$Gene, levels=c('GTC',tophits))
d$group <- factor(substr(d$Sample,0,1), levels=c('T','D','P'))
d$Sample <- d$Sample %>% stringr::str_replace(., 'T3', 'T1') %>% stringr::str_replace(., 'T5', 'T2') %>% stringr::str_replace(., 'T6', 'T3')
d$Sample <- factor(d$Sample, levels=d$Sample[order(d$group, decreasing=T)] %>% unique())
sg.order <- sgRes[which(sgRes$sgRNA %in% top4guide),]
sg.order$Gene <- factor(sg.order$Gene, levels=tophits)
sg.order <- sg.order$sgRNA[with(sg.order, order(Gene, -logFC))]
d$sgRNA <- factor(d$sgRNA, levels=c('GTC', sg.order))

# set colors for each sgRNA
col <- structure(c('gray40',ggcolor(6)), names=levels(d$Gene))

# ensure same ylimits are used for all plots
ylim <- c(floor(range(d$Value, na.rm=T)[1]), ceiling(range(d$Value, na.rm=T)[2]))
p <- lapply(levels(d$Gene)[-1], function(g){
    ifelse(g=='PTEN_ZBTB7B', y.axis <- 'sgRNA z-score', y.axis <- '')
    ggplot(d[which(d$Gene %in% c(g,'GTC')),], aes(x=Sample, y=Value, shape=sgRNA, color=Gene, group=sgRNA)) + 
        geom_point(size=0.6) + theme_test() + # geom_line(alpha=0.8, size=0.3)  +
        ylim(ylim) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=6), legend.position=c(0.12,.75)) +
        scale_color_manual(values=col, guide='none') + addSmallLegend(1.2) + guides(color = "none") +
        labs(title=g, x='', y=y.axis) 
})
p <- cowplot::plot_grid(plotlist=p, align='hv', ncol=3)
ggsave(plot=p, 'Figs/plots_tophits_by_sample_byGene_v2.1.pdf', width=9, height=3)

