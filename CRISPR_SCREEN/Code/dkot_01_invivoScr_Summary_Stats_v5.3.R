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

ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
options("ggrastr.default.dpi" = 750)


source('Code/samba_v2.0.R')
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
rc <- merge(read.delim('Data/dkot_invivoScrTumor_v2m2.txt'),
            read.delim('Data/dkot_invivoScrSpleen_v2m2.txt')[,c(1,3:13)],
            by='sgRNA')

colnames(rc) <- colnames(rc) %>% stringr::str_replace(., 'T3', 'T1') %>% stringr::str_replace(., 'T5', 'T2') %>% stringr::str_replace(., 'T6', 'T3') %>% 
    stringr::str_replace(., 'S3', 'S1') %>% stringr::str_replace(., 'S5', 'S2') %>% stringr::str_replace(., 'S6', 'S3') %>% 
    stringr::str_replace(., 'D3', 'D1') %>% stringr::str_replace(., 'D5', 'D2') %>% stringr::str_replace(., 'D6', 'D3')
rc[is.na(rc)] <- 0
libdata <- read.delim('Data/dkot_libaryInformation.txt')
sampledata <- data.frame(samples = colnames(rc[,3:ncol(rc)]),
                         Donor = substr(colnames(rc[,3:ncol(rc)]),2,2),
                         Tissue = c(S = 'Spleen', P = 'Plasmid', D = 'Culture',T = 'Tumor')[substr(colnames(rc[,3:ncol(rc)]),1,1)],
                         Type = c('Screen','Cell','Plasmid')[grepl('Cell',colnames(rc[,3:ncol(rc)])) + 
                                                             (2*grepl('Plasmid',colnames(rc[,3:ncol(rc)]))) + 1])
sampledata$label <- c(paste0('D',sampledata$Donor[1:11],' ', 'Tumor', substr(sampledata$sample[1:11],3,3)), paste0('D',1:3, ' Cell'), 
                      'Plasmid', c(paste0('D',sampledata$Donor[16:26],' ', 'Spleen', substr(sampledata$sample[16:26],3,3))))


## Show library distribution by gRNA type using ridgeplots.
df <-  cbind(rc[,1:2], apply(rc[,3:ncol(rc)],2,function(x) 1e6*x/sum(x)))
df$Gene[!(df$Gene %in% c('NTC','GTC'))] <- 'gRNA'
df <- reshape2::melt(df[,-1])
df$Sample <- as.character(df$variable)
df <- merge(df,sampledata, by.y='samples',by.x='Sample',all.x=T,all.y=F)
df$Tissue[grep('Cell', df$Sample)] <- 'Cell'
df$Tissue <- factor(df$Tissue, levels = c('Plasmid','Culture','Tumor','Spleen','Cell'))
df$value <- log2(1+df$value)
p1 <- ggplot(df, aes(x = value, fill = Gene, y = Tissue)) + ggridges::geom_density_ridges2(alpha=0.5) + 
    xlim(c(-1,7.5)) + cowplot::theme_cowplot() + addSmallLegend(4) +
    labs(x = 'gRNA density (log2 cpm)', y = 'Sample type')
p2 <- ggplot(df, aes(x = value, y = Gene, fill = Tissue)) + ggridges::geom_density_ridges2(alpha=0.5) + 
    xlim(c(-1,7.5)) + cowplot::theme_cowplot() + addSmallLegend(4) +
    labs(x = 'gRNA density (log2 cpm)', y = 'Guide type')
p <- p1 + p2 + plot_annotation('gRNA density by guide and sample type',theme=theme(plot.title=element_text(hjust=0.5)))
ggsave(plot = p, filename = file.path(paste0('Figs/qc_invivo_countDensity_plot_v3.pdf')), height = 4, width = 8)


## Sample correlation heatmaps
lcpm <- log2(1 + apply(rc[,3:ncol(rc)],2,function(x) 1e6*x/sum(x)))
corr <- cor(lcpm, method = 's')
dimnames(corr) <- list(sampledata$label, sampledata$label)
pheatmap::pheatmap(corr, scale = 'none', treeheight_row = 12, cluster_cols = T, 
                   treeheight_col = 12,  clustering_method = 'ward.D2', 
                   angle_col = 90, width = 4.7, height = 4.5, 
                   cutree_rows = 1, cutree_cols = 1,
                   filename = 'Figs/qc_invivo_corr_heatmap_v3.pdf')
write.table(corr, file = 'Data/qc_invivo_corr_heatmap_v3.txt',sep = '\t', quote = F)


## Read count distribution
lcpm <- log2(1 + apply(rc[,3:ncol(rc)],2,function(x) 1e6*x/sum(x)))
df <- reshape2::melt(lcpm)
colnames(df) <- c('guide','sample','value')
df <- merge(df,sampledata, by.y='samples',by.x='sample',all.x=T,all.y=F)
df.summary <- df %>% group_by(sample) %>% reframe(Donor = unique(Donor), value = mean(value))
df.summary$Donor <- factor(df.summary$Donor, levels=c('l','1','2','3'))
df.summary$Tissue <- factor(substr(df.summary$sample,1,1), levels=c('P','D','S','T'))
#df.summary <- with(df.summary, order(Donor,-value, decreasing=F))
y.order <- df.summary$sample[with(df.summary, order(Donor, Tissue, -value, decreasing=F))]
df$sample <- factor(df$sample, levels = rev(y.order))
df$Donor <- factor(df$Donor, levels=c('l','1','2','3'))
pal <- structure(ggcolor(4), names = c('l','1','2','3'))
p <- ggplot(df, aes(x = sample, y = value, color = Donor)) + 
    ggrastr::geom_jitter_rast(size = 0.1, pch = 19, alpha= 0.3) + 
    geom_boxplot(color = 'gray50',outlier.size = 0.2,outlier.alpha= 0.7)  +
    coord_flip() + labs(x = '', y = 'sgRNA counts (log2 cpm)', title = 'Library distribution') + 
    cowplot::theme_cowplot() + addSmallLegend(4)
ggsave(plot = p, filename = file.path(paste0('Figs/qc_invivo_countDist_v3.png')), height = 6, width = 4)



## Cumulative prob functions for all groups
y.order <- sampledata[with(sampledata, order(factor(sampledata$Tissue, levels = c('Plasmid','Culture','Spleen','Tumor')), samples)),]
cdf <- lapply(y.order$samples, function(s) { ecdf(df[which(df$sample == s),'value']) })# %>% structure(., names = unique(d$group))
for(i in 1:length(cdf)) p1 <- p1 + geom_function(aes(color = y.order$Tissue[i]), fun = cdf[[i]])
cdfs <- lapply(1:length(cdf), function(i) geom_function(aes(color = y.order$Tissue[i]), fun = cdf[[i]]))
p1 <- ggplot() + xlim(range(df$value)) + cdfs + labs(y = 'Cumulative probability', x = 'sgRNA counts (log2 cpm)') + cowplot::theme_cowplot() +
    theme(legend.position = c(0.6,0.25)) + guides(color=guide_legend(title='Samples')) 
cdf <- lapply(unique(y.order$Tissue), function(group) { ecdf(df[which(df$Tissue == group),'value']) })# %>% structure(., names = unique(d$group))
p2 <- ggplot() + xlim(range(df$value)) + 
    geom_function(aes(color = 'Plasmid'), fun = cdf[[1]]) +
    geom_function(aes(color = 'Culture'), fun = cdf[[2]]) +
    geom_function(aes(color = 'Spleen'), fun = cdf[[3]]) +
    geom_function(aes(color = 'Tumor'), fun = cdf[[4]]) +
    labs(y = 'Cumulative probability', x = 'sgRNA counts (log2 cpm)') + cowplot::theme_cowplot() +
    theme(legend.position = c(0.6,0.25)) + guides(color=guide_legend(title='Sample type')) 
p <- p1 + p2 + plot_annotation('ECDF plots for gRNA library distribution',theme=theme(plot.title=element_text(hjust=0.5)))
ggsave(plot = p, filename = file.path(paste0('Figs/qc_invivo_cdf_v3.pdf')), height = 4, width = 6)


## KS test for differences in library distribution 
dge <- readRDS(file = paste0('Data/dkot_invivoScr_dge.rds'))
lcpm <- log2(1 + apply(dge$counts,2,function(x) 1e6*x/sum(x)))
lcpm <- dge$counts
lcpm <- data.frame(Plasmid=lcpm[,15], Culture=rowSums(lcpm[,12:14]),Tumor=rowSums(lcpm[,1:11]))
d <- lcpm
ks <- lapply(c('Plasmid','Culture','Tumor'), function(i) { 
    lapply(c('Plasmid','Culture','Tumor'), function(j) { 
        ks.test(x=d[,i], y=d[,j],alternative='t', simulate.p.value=T, B=1000)$p.value
        }) %>% unlist()
    })%>% unlist()
ks <- matrix(ks, byrow=T, nrow=3, dimnames=list(c('Plasmid','Culture','Tumor'),c('Plasmid','Culture','Tumor')))



## Calculate Gini index for all samples
det <- edgeR::gini(lcpm)
det <- data.frame(sample = names(det), 
                  Tissue = sampledata$Tissue[match(names(det),sampledata$samples)],
                  Gini = det)
det$Tissue <- factor(det$Tissue, levels = c('Plasmid','Culture','Cell','Spleen','Tumor'))
y.order <- det$sample[with(det, order(Tissue,-Gini, decreasing = F))]
det$sample <- factor((det$sample), levels = rev(y.order))
p1 <- ggplot(det, aes(x = sample, y = Gini, fill = Tissue)) + 
    geom_bar(stat = 'identity') + coord_flip() + theme_classic() + 
    labs(x = '', y = 'Gini coefficient', title = 'Readout diversity') + theme(legend.position = 'none') 
d <- lcpm
d[d > 0] <- 1
det <- edgeR::gini(d)
det <- data.frame(sample = names(det), 
                  Tissue = sampledata$Tissue[match(names(det),sampledata$samples)],
                  Gini = det)
det$Tissue <- factor(det$Tissue, levels = c('Plasmid','Culture','Cell','Spleen','Tumor'))
det$sample <- factor((det$sample), levels = rev(y.order))
p2 <- ggplot(det, aes(x = sample, y = Gini, fill = Tissue)) + 
    geom_bar(stat = 'identity') + coord_flip() + theme_classic() + 
    labs(x = '', y = 'Gini coefficient', title = 'ZeroCount diversity') + theme(legend.position = 'none') 
ggsave(plot = p1 | p2, filename = 'Figs/qc_invivo_gini_v3.pdf', height = 4, width = 5)




## sgRNA-detection rates
d <- 100*apply(rc[,3:ncol(rc)], 2, function(x) sum(x > 0))/nrow(rc)
det <- data.frame(sample = names(d), 
                  Tissue = sampledata$Tissue[match(names(d), sampledata$samples)],
                  Donor = substr(names(d),2,2) %>% factor(levels = c('l','1','2','3')),
                  detection = d)
det$Tissue <- factor(det$Tissue, levels=c('Plasmid','Culture','Spleen','Tumor'))
det$label <- c(paste0('D',sampledata$Donor[1:11],' ', 'Tumor', substr(sampledata$sample[1:11],3,3)), paste0('D',1:3, ' Cell'), 
                      'Plasmid', c(paste0('D',sampledata$Donor[16:26],' ', 'Spleen', substr(sampledata$sample[16:26],3,3))))
y.order <- det$sample[with(det, order(Donor, Tissue, -detection, decreasing = F))]
det$sample <- factor((det$sample), levels = rev(y.order))
p <- ggplot(det, aes(x = sample, y = detection, fill = Donor)) + 
    geom_bar(stat = 'identity') + coord_flip() + theme_classic() +
    scale_x_discrete(labels=det$label[rev(with(det, order(Donor, Tissue, -detection, decreasing = F)))]) + 
    labs(x = '', y = 'Library detection (%)', title = 'gRNA-detection rates') + theme(legend.position = 'none')
ggsave(plot = p, filename = 'Figs/qc_invivo_DetRate_v3.pdf', height = 4, width = 3)








