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
SelectLabels <- function(x, fc.threshold = 1, keeplist=''){
    x <- x[which((x$FDR < 0.05) & (abs(x$logFC) > fc.threshold)),]
    x <- x[order(-x$F),]
    topsig.up <- x$transcript_id[which(x$logFC > 0)] %>% head(n=8)
    topsig.dn <- x$transcript_id[which(x$logFC < 0)] %>% head(n=8)
    x <- x[order(-abs(x$logFC)),]
    topfc.up <- x$transcript_id[which(x$logFC > 0)] %>% head(n=8)
    topfc.dn <- x$transcript_id[which(x$logFC < 0)] %>% head(n=8)
    genes <- unique(c(topsig.up, topsig.dn, topfc.up, topfc.dn, keeplist))
    genes <- genes[which(!is.null(genes) & genes != '')]
    return(genes)
}


## load prepped data
data <- read.delim('Data/02a_dkotrna_txcounts.txt', row.names=1)
meta <- read.delim('Data/02b_dkotrna_metadata_table.txt')
rownames(meta) <- meta$sample
meta$sample <- factor(meta$sample, levels=meta$sample)
meta$genotype <- factor(meta$genotype, levels=c('GTC','NR4A1','SOCS3','DKO'))
meta$time <- factor(meta$time, levels=c('0','6','24','72'))
gtf <- read.delim('Reference/gencode_GRCh38.96_tx_info.txt')
tx2use <- gtf$transcript_id[which(gtf$transcript_support_level < 4)]

###################

## get lcpm data from all samples
meta$groups <- paste0(meta$genotype,'_', meta$time)
dge <- DGEList(counts = data[,meta$sample], genes=gtf, group = factor(meta$groups, levels=unique(meta$groups)), remove.zeros = T)
dge <- dge[which(rownames(dge) %in% tx2use),]
## Pre-process
dge <- dge[filterByExpr(dge, group=dge$samples$group, min.count = 5, min.total.count = 10),, keep.lib.sizes=F]
dge <- calcNormFactors(object = dge, method = "TMMwsp")
design <- model.matrix(~ 0 + groups, data=meta)
colnames(design) <- stringr::str_remove(colnames(design), 'groups') 
dge <- estimateGLMCommonDisp(dge, design = design)
dge <- estimateGLMTrendedDisp(dge, method = 'bin.loess', span = 1/3, design = design)
dge <- estimateGLMTagwiseDisp(dge, design = design)
# Fit model and run QLF test
fit <- glmQLFit(dge, design = design, robust=T)
lcpm <- data.frame(gtf[rownames(fit$coefficients), c(3,2)], edgeR::cpm(dge, log=T))
write.table(lcpm, 'Data/03a_lcpm/03a_lcpm_tx.txt', sep='\t', quote=F)


###################

## Run pairwise analysis of genotypes, separately for each timepoint
for(i in 1:3){
    for(j in (i+1):4){
        genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
        genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
        for(timepoint in c(0,6,24,72)){
            samples.to.use <- meta$sample[which((meta$time == timepoint) & (meta$genotype%in%c(genotype1,genotype2)))]
            dge <- DGEList(counts = data[,samples.to.use], genes=gtf, remove.zeros = T, group = factor(meta[samples.to.use, 'genotype']))
            dge <- dge[which(rownames(dge) %in% tx2use),]
            
            ## Pre-process
            #filter low-expressed transcripts
            group <- dge$samples$group
            lcpm.raw <- cpm(dge, log = T)
            dge <- dge[filterByExpr(dge, group=dge$samples$group, min.count = 5, min.total.count = 10),, keep.lib.sizes=F]
            lcpm.filt <- cpm(dge, log = T)
            #get cutoff
            L <- mean(dge$samples$lib.size) * 1e-6
            M <- median(dge$samples$lib.size) * 1e-6
            lcpm.cutoff <- log2(2/M + 1/L)
            #plot filtered vs unfiltered
            p1 <- ggplot(reshape2::melt(lcpm.raw)) + geom_density(aes(x = value, color = Var2)) + 
                geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
                labs(title = 'Raw data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
                theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
                      legend.title = element_text(size=8),legend.text = element_text(size=6), 
                      legend.background = element_rect(fill = "white", color = "black")) 
            p2 <- ggplot(reshape2::melt(lcpm.filt)) + geom_density(aes(x = value, color = Var2)) + 
                geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
                labs(title = 'Filtered data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
                theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
                      legend.title = element_text(size=8),legend.text = element_text(size=6), 
                      legend.background = element_rect(fill = "white", color = "black")) 
            rm(lcpm.cutoff, L,M, lcpm.raw, lcpm.filt)
            
            
            ## Normalize
            dge.unnorm <- cpm(dge, log = T, normalized.lib.sizes = F)
            dge <- calcNormFactors(object = dge, method = "TMMwsp")
            dge.norm <- cpm(dge, log = T, normalized.lib.sizes = T)
            #visualize distributions with Norm cpm
            p3 <- ggplot(reshape2::melt(dge.unnorm)) + geom_boxplot_jitter(aes(y = value,x = Var2, color = Var2), size = 0.2) + 
                labs(title = 'Unnormalized data',y = 'Log-cpm', x = '', color = 'Sample') + theme_classic() + 
                theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
            p4 <- ggplot(reshape2::melt(dge.norm)) + geom_boxplot_jitter(aes(y = value,x = Var2, color = Var2), size = 0.2) + 
                labs(title = 'Normalized data',y = 'Log-cpm', x = '', color = 'Sample') + theme_classic() + 
                theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
            ggsave(plot = (p1 + p2) / (p3 + p4), filename = paste0('Figures/03a_qc_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.pdf'), height = 6, width = 5)
            
            
            # Estimate dispersion
            geno <- as.integer(dge$samples$group==levels(meta$genotype)[max(as.integer(meta[samples.to.use,]$genotype))])
            design <- model.matrix(~ geno)
            dge <- estimateGLMCommonDisp(dge, design = design)
            dge <- estimateGLMTrendedDisp(dge, method = 'bin.loess', span = 1/3, design = design)
            dge <- estimateGLMTagwiseDisp(dge, design = design)
            
            
            # Fit model and run QLF test
            fit <- glmQLFit(dge, design = design, robust=T)
            lcpm <- data.frame(gtf[rownames(fit$coefficients), c(3,2)], edgeR::cpm(dge, log=T))
            write.table(lcpm, paste0('Data/03a_lcpm/03a_lcpm_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'), sep = '\t', col.names = T, quote = F)
            de <- topTags(glmQLFTest(fit, coef = 'geno'), n = Inf)$table
            de <- de[order(de[,5], decreasing=T),]
            write.table(de, paste0('Data/03b_de/03b_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
        }
    }
}


## make volcanos 
for(i in 1:3){
    for(j in (i+1):4){
        genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
        genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
        p <- lapply(c(0,6,24,72), function(timepoint){
                de <- read.delim(paste0('Data/03b_de/03b_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))

                ## make volcano plots
                genes2label <- SelectLabels(de)
                p <- ggplot(de, aes(x=logFC, y=-log10(FDR))) +
                    ggrastr::geom_point_rast(color='gray20', size=0.7) +
                    ggrastr::geom_point_rast(color='gray65', size=0.2) +
                    ggrastr::geom_point_rast(data=de[which(de$logFC < -1 & de$FDR < 0.05),], color='steelblue', size=0.4) +
                    ggrastr::geom_point_rast(data=de[which(de$logFC >  1 & de$FDR < 0.05),], color='firebrick', size=0.4) +
                    ggrepel::geom_text_repel(data=de[which(de$transcript_id %in% genes2label),], aes(label=transcript_name), size=1.8, min.segment.length=0, nudge_y=0.5, segment.size=0.2, segment.alpha=0.5) + 
                    labs(x='Log2 fold-change', y='q value (-log10)', title=paste0('DE ',timepoint,'hrs')) + 
                    theme_test()
                return(p)
            })
            ggsave(plot = cowplot::plot_grid(plotlist=p, nrow=1, align='hv'), 
                   filename = paste0('Figures/03b_volcano_tx_',genotype1,'-',genotype2,'_v1.3.pdf'), height = 2.5, width = 11)
    }
}

###################


## Run epistasis analysis, separately for each timepoint
dm <- data.frame(GTC=as.integer(meta$genotype %in% c('GTC')),
                 NR4A1=as.integer(meta$genotype %in% c('NR4A1','DKO')),
                 SOCS3=as.integer(meta$genotype %in% c('SOCS3','DKO')),
                 DKO=as.integer(meta$genotype=='DKO'),row.names=meta$sample)
for(timepoint in as.character(c(0,6,24,72))){
            samples.to.use <- meta$sample[which(meta$time == timepoint)]
            dge <- DGEList(counts = data[,samples.to.use], genes=gtf, 
                           group = factor(meta[samples.to.use, 'genotype']), remove.zeros = T)
            dge <- dge[which(rownames(dge) %in% tx2use),]
            dge$dm <- dm[samples.to.use, ]
            
            ## Pre-process
            #filter low-expressed transcripts
            group <- dge$samples$group
            lcpm.raw <- cpm(dge, log = T)
            dge <- dge[filterByExpr(dge, group=dge$samples$group),, keep.lib.sizes=F]
            lcpm.filt <- cpm(dge, log = T)
            #get cutoff
            L <- mean(dge$samples$lib.size) * 1e-6
            M <- median(dge$samples$lib.size) * 1e-6
            lcpm.cutoff <- log2(2/M + 1/L)
            #plot filtered vs unfiltered
            p1 <- ggplot(reshape2::melt(lcpm.raw)) + geom_density(aes(x = value, color = Var2)) + 
                geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
                labs(title = 'Raw data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
                theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
                      legend.title = element_text(size=8),legend.text = element_text(size=6), 
                      legend.background = element_rect(fill = "white", color = "black")) 
            p2 <- ggplot(reshape2::melt(lcpm.filt)) + geom_density(aes(x = value, color = Var2)) + 
                geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
                labs(title = 'Filtered data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
                theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
                      legend.title = element_text(size=8),legend.text = element_text(size=6), 
                      legend.background = element_rect(fill = "white", color = "black")) 
            rm(lcpm.cutoff, L,M, lcpm.raw, lcpm.filt)
            
            
            ## Normalize
            dge.unnorm <- cpm(dge, log = T, normalized.lib.sizes = F)
            dge <- calcNormFactors(object = dge, method = "TMMwsp")
            dge.norm <- cpm(dge, log = T, normalized.lib.sizes = T)
            #visualize distributions with Norm cpm
            p3 <- ggplot(reshape2::melt(dge.unnorm)) + 
                geom_boxplot_jitter(aes(y = value,x = Var2, color = Var2), size = 0.2) + 
                labs(title = 'Unnormalized data',y = 'Log-cpm', x = '', color = 'Sample') + theme_classic() + 
                theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
            p4 <- ggplot(reshape2::melt(dge.norm)) + 
                geom_boxplot_jitter(aes(y = value,x = Var2, color = Var2), size = 0.2) + 
                labs(title = 'Normalized data',y = 'Log-cpm', x = '', color = 'Sample') + theme_classic() + 
                theme(legend.position = 'none',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
            ggsave(plot = (p1 + p2) / (p3 + p4), 
                   filename = paste0('Figures/03c_qc_tx_',timepoint,'hrs_','epistasis','.pdf'), height = 6, width = 5)
            
            
            # Estimate dispersion
            design <- model.matrix(~ NR4A1 + SOCS3 + DKO, data=dge$dm)
            dge <- estimateGLMCommonDisp(dge, design = design)
            dge <- estimateGLMTrendedDisp(dge, method = 'bin.loess', span = 1/3, design = design)
            dge <- estimateGLMTagwiseDisp(dge, design = design)
            
            
            # Fit model and run QLF test
            fit <- glmQLFit(dge, design = design)
            lcpm <- data.frame(gtf[rownames(fit$coefficients), c(3,2)], edgeR::cpm(dge, log=T))
            write.table(lcpm, paste0('Data/03d_lcpm_epistasis/03c_lcpm_tx_',timepoint,'hrs_epistasis.txt'), 
                        sep = '\t', col.names = T, quote = F)
            de <- topTags(glmQLFTest(fit, coef='DKO'), n = Inf)$table
            write.table(de, paste0('Data/03d_de_epistasis/03d_de_tx_',timepoint,'hrs_epistasis.txt'), 
                        sep = '\t', col.names = T, row.names = F, quote = F)

}


## make volcano plots
p <- lapply(as.character(c(0,6,24,72)), function(timepoint){
    de <- read.delim(paste0('Data/03d_de_epistasis/03d_de_tx_',timepoint,'hrs_epistasis.txt'))
    genes2label <- SelectLabels(de, 1)
    p <- ggplot(de, aes(x=logFC, y=-log10(FDR))) +
        ggrastr::geom_point_rast(color='gray20', size=1) +
        ggrastr::geom_point_rast(color='gray65', size=0.6) +
        ggrastr::geom_point_rast(data=de[which(de$logFC < -1 & de$FDR < 0.05),], color='steelblue', size=0.4) +
        ggrastr::geom_point_rast(data=de[which(de$logFC >  1 & de$FDR < 0.05),], color='firebrick', size=0.4) +
        ggrepel::geom_text_repel(data=de[which(de$transcript_id %in% genes2label),], aes(label=transcript_name), size=2.5, min.segment.length=0, nudge_y=0.5, segment.size=0.2, max.overlaps=40, segment.alpha=0.5) + 
        labs(x='Log2 fold-change', y='q value (-log10)', title=paste0('DE ',timepoint,'hrs epistasis')) +
        theme_test()
    return(p)
})
ggsave(plot = cowplot::plot_grid(plotlist=p, nrow=1, align='hv'), height = 3, width = 15, 
       filename = paste0('Figures/03d_volcano_tx_epistasis_v1.3.pdf'))


