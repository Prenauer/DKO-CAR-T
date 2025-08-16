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


## Accessory functions 
SelectLabels <- function(x, fc.threshold = 1, keeplist=''){
    x <- x[which((x$FDR < 0.05) & (abs(x$logFC) > fc.threshold)),]
    x <- x[order(x$PValue),]
    topsig.up <- x$transcript_id[which(x$logFC > 0)] %>% head(n=8)
    topsig.dn <- x$transcript_id[which(x$logFC < 0)] %>% head(n=8)
    x <- x[order(-abs(x$logFC)),]
    topfc.up <- x$transcript_id[which(x$logFC > 0)] %>% head(n=8)
    topfc.dn <- x$transcript_id[which(x$logFC < 0)] %>% head(n=8)
    genes <- unique(c(topsig.up, topsig.dn, topfc.up, topfc.dn, keeplist))
    genes <- genes[which(!is.null(genes) & genes != '')]
    return(genes)
}

Optimal_k <- function(x, k.range = c(3:10)) {
    d <- dist(t(scale(t(x))))
    hc <- hclust(d, method='ward.D2')
    asw <- sapply((min(k.range)-1):(max(k.range)+1), function(k) mean(cluster::silhouette(cutree(hc, k), d)[,3]) )
    plot(1:length(asw),asw)
    asw.diff <- sapply(2:(length(asw)-1), function(i) asw[i]-mean(asw[(i-1):(i+1)]))
    plot(k.range,asw.diff)
    optim.k <- k.range[which(asw.diff==max(asw.diff))]
    return(optim.k)
}

GroupMeans <- function(x){
    y=reframe(data.frame(meta[colnames(x),2:3], scale(t(x))), .by=c('genotype','time'), across(where(is.numeric), mean))
    cnames <- paste0(y[,1], '_', y[,2])
    y <- t(y[,-c(1:2)]) %>% data.frame()
    colnames(y) <- cnames
    return(y)
}
DetectNonRedundantTerms <- function(x){
    x <- stringr::str_split(x, ',')
    tmp <- do.call(rbind, lapply(1:(length(x) - 1), function(i){
        tmp <- lapply((i+1):length(x), function(j) sum(x[[i]] %in% x[[j]])/length(x[[j]]))
        return(data.frame(pw1=i, pw2=(i+1):length(x), sim=unlist(tmp)))
    }))
    tmp <- do.call(rbind, lapply(1:length(x), function(i){
        tmp <- lapply(c(1:length(x))[-i], function(j) sum(x[[i]] %in% x[[j]])/length(x[[j]]))
        return(data.frame(pw1=i, pw2=c(1:length(x))[-i], sim=unlist(tmp)))
    }))
    
    rm.list <- tmp$pw2[which(tmp$sim == 1)] %>% unique()
    return(!(c(1:length(x)) %in% rm.list))
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


### Full heatmap of DEGs
# get de genes from all comparisons
deg <- (lapply(1:3, function(i){
    deg <- (lapply((i+1):4, function(j){
        genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
        genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
        deg <- lapply(c(0,6,24,72), function(timepoint){
            de <- read.delim(paste0('Data/03b_de/03b_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))
            return(de$transcript_id[which(de$FDR < 0.01 & abs(de$logFC) > 0.5)])
        }) %>% unlist() %>% unique()
        return(deg)
    }))
    return(deg)
})) 
deg <- unique(unlist(deg))  
# Make heatmap of lcpm values
lcpm <- read.delim('Data/03a_lcpm/03a_lcpm_tx.txt')
k <- Optimal_k(lcpm[which(lcpm$transcript_id %in% deg), 3:ncol(lcpm)])
hm.filename <- paste0('Figures/04a_hm_deg_allConditions.pdf')
hm <- pheatmap::pheatmap(lcpm[deg, 3:ncol(lcpm)], scale='row', cluster_cols=T, labels_row = gtf[deg, 'transcript_name'],  
                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
                         show_rownames=F, treeheight_col=15, treeheight_row=20, clustering_method='ward.D2', annotation_col = meta[,3, drop=F], 
                         angle_col='90', 
                         cutree_rows=k, cutree_cols=5, filename=hm.filename, height=7, width=7)


####################
pw <- list()
pw$hallmark <- fgsea::gmtPathways('Reference/h.all.v2023.1.Hs.symbols.gmt.txt')
pws2exclude <- c('HALLMARK_COAGULATION','HALLMARK_MYOGENESIS','HALLMARK_SPERMATOGENESIS','HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_BILE_ACID_METABOLISM','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_PANCREAS_BETA_CELLS','HALLMARK_XENOBIOTIC_METABOLISM','HALLMARK_HEME_METABOLISM','HALLMARK_UV_RESPONSE_UP','HALLMARK_UV_RESPONSE_DN')
pw$hallmark <- pw$hallmark[!(names(pw$hallmark) %in% pws2exclude)]
pw$tf <- fgsea::gmtPathways('Reference/c3.TF_targets.v7.4.symbols.gmt.txt')

### Analyze by timepoint (Hallmark)
timepoints <- c(0,6,24,72)
for(timepoint in timepoints){
    # get de genes from all comparisons
    deg <- (lapply(1:3, function(i){
        deg <- (lapply((i+1):4, function(j){
            genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
            genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
            de <- read.delim(paste0('Data/03b_de/03b_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))
            deg <- de$transcript_id[which(de$FDR < 0.01 & abs(de$logFC) > 0.5)]
            return(deg)
        }))
        return(deg)
    })) 
    deg <- unique(unlist(deg))  
    
    # get log-norm data for comparisons
    lcpm <- read.delim(paste0('Data/03a_lcpm/03a_lcpm_tx_',timepoint,'hrs.txt'))

    # Make heatmap
    k <- Optimal_k(lcpm[which(lcpm$transcript_id %in% deg), 3:ncol(lcpm)], 3:5)
    windsorize <- function(x) {
        ul <- 2 
        ll <- -2 
        x[which(x > ul)] <- ul
        x[which(x < ll)] <- ll
        return(x)
    } 
    d <- t(scale(t(lcpm[deg, 3:ncol(lcpm)])))
    d <- apply(d, 2, function(x) windsorize(x))
    hm <- pheatmap::pheatmap(d, scale='none', cluster_cols=T, labels_row = gtf[deg, 'transcript_name'],  
                             color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),
                             show_rownames=F, treeheight_col=15, treeheight_row=20, clustering_method='ward.D2', annotation_legend=F,
                             angle_col='90', main=paste0(timepoint, 'hr DEG'), annotation_col = meta[which(meta$time==timepoint), 2, drop=F],
                             cutree_rows=k, filename=NA, height=4, width=4)
    # collect and save heatmap clustering info
    cl <- cutree(hm$tree_row, k)[hm$tree_row$order]
    cl <- structure(match(1:k, unique(cl))[cl], names=names(cl)) # renumber clusters
    cl <- merge(gtf[,1:3],data.frame(transcript_id=names(cl), clust=cl), by='transcript_id')[,3:4] %>% unique()
    write.table(cl, paste0('Data/04a_de_clusters/04a_clust_deg_timepoint_',timepoint,'hrs_hallmark.txt'), sep='\t', quote=F, row.names=F)
    
    
    
    ### Cluster Pathway analyses  
    universe <- unique(gtf$gene_name[which(!is.na(gtf$gene_name))])
    # run ORA and save
    res <- do.call(rbind, lapply(names(pw), function(n){
        do.call(rbind, lapply(unique(cl$clust), function(clust){
            res <- fgsea::fora(pw[[n]], genes = cl$gene_name[which(cl$clust==clust)], universe = universe, minSize=2, maxSize=1000)
            if(nrow(res) > 0) res <- data.frame(pathway_set=n, clust=clust, n_clust_genes=sum(cl$clust==clust), res)
            return(res)
        }))
    }))
    res$pathway <- stringr::str_remove(res$pathway, 'REACTOME_') %>% stringr::str_remove(., 'HALLMARK_') %>% stringr::str_remove(., 'GOBP_') %>% stringr::str_replace_all(., '_', ' ') 
    res <- res[which(res$overlap > 1),]
    res <- res[with(res, order(pathway_set, clust, pval)),]
    res$recall <- res$overlap / res$size
    res$precision <- res$overlap / res$n_clust_genes
    res$overlapGenes <- sapply(res$overlapGenes, function(x) paste0(x, collapse=','))
    write.table(res, paste0('Data/04b_pathway_analysis/04b_gm_ora_timepoint_',timepoint,'hrs.txt'), sep='\t', quote=F, row.names=F)
    
    ### Prepare pw bubble plots
    # make gene lists
    res <- res[which((res$overlap > 4) & (res$padj < 0.05) & (res$size <= 600)),]
    res$NonredundantTerms <- T
    res$NonredundantTerms[which(res$pathway_set=='bp')] <- DetectNonRedundantTerms(res$overlapGenes[which(res$pathway_set=='bp')])
    res <- res[which(res$NonredundantTerms),]
    res$overlapGenes <- strsplit(res$overlapGenes, ',')
    p.data <- slice_head(res, by=c('pathway_set','clust'), n=4)
    p.data$label <- paste0(toupper(substr(p.data$pathway,0,1)),  substr(tolower(p.data$pathway),2,80)) %>% stringr::str_wrap(., 30) #, ' (q=', formatC(res$padj, 2, format='e'),')') 
    p.data$pw_ordered <- paste0(rev(formatC(01:nrow(p.data), width=2, flag='0')), '_', p.data$label) %>% factor()

    p1 <- ggdraw() + draw_grob(hm$gtable)
    p2 <- ggplot(p.data[which(p.data$pathway_set=='bp'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
        geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='bp')]), position='right') + theme_light() +
        labs(title='BP', x= 'q value (-log10)', y='', size='Recall', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
    p3 <- ggplot(p.data[which(p.data$pathway_set=='hallmark'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
        geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='hallmark')]), position='right') + theme_light() +
        labs(title='Hallmark', x= 'q value (-log10)', y='', size='Recall', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
    p4 <- ggplot(p.data[which(p.data$pathway_set=='tf'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
        geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='tf')]), position='right') + theme_light() +
        labs(title='TF', x= 'q value (-log10)', y='', size='Recall', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
    p <- plot_grid(p2, p3, p4, align='v', rel_widths = c(1,1,1), axis='none', ncol=3)
    p <- plot_grid(p1, p, align='v', rel_widths = c(0.2,1), axis='b', ncol=2)
    ggsave(plot=p, paste0('Figures/04b_bubl_pw_timepoint_',timepoint,'hrs_hallmark.pdf'), height=6, width=15)
}


## Make new heatmaps
timepoints <- c(0,6,24,72)
p <- lapply(timepoints, function(timepoint){
  # get de genes from all comparisons
  deg <- (lapply(1:3, function(i){
    deg <- (lapply((i+1):4, function(j){
      genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
      genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
      de <- read.delim(paste0('Data/03b_de/03b_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))
      deg <- de$transcript_id[which(de$FDR < 0.01 & abs(de$logFC) > 0.5)]
      return(deg)
    }))
    return(deg)
  })) 
  deg <- unique(unlist(deg))  
  
  # get log-norm data for comparisons
  lcpm <- read.delim(paste0('Data/03a_lcpm/03a_lcpm_tx_',timepoint,'hrs.txt'))
  
  k <- Optimal_k(lcpm[which(lcpm$transcript_id %in% deg), 3:ncol(lcpm)], 3:5)
  windsorize <- function(x) {
    ul <- 2 
    ll <- -2 
    x[which(x > ul)] <- ul
    x[which(x < ll)] <- ll
    return(x)
  } 
  d <- t(scale(t(lcpm[deg, 3:ncol(lcpm)])))
  d <- apply(d, 2, function(x) windsorize(x))
  hm <- pheatmap::pheatmap(d, scale='none', cluster_cols=F, labels_row = gtf[deg, 'transcript_name'],  show_colnames=F,
                           show_rownames=F, clustering_method='ward.D2', annotation_legend=F,
                           cutree_rows=k, filename=NA)
  # create row annotation for heatmap
  cl <- cutree(hm$tree_row, k)[hm$tree_row$order]
  cl <- structure(match(1:k, unique(cl))[cl], names=names(cl)) # renumber clusters
  row.anno <- data.frame(C=cl[rownames(d)], row.names=rownames(d))
  row.anno$C <- paste0(formatC(timepoint, width=2, flag='0'),'h-',row.anno$C)
  color.anno <- list(C=structure(c(ggcolor(4),'orange')[1:length(unique(row.anno$C))], names=sort(unique(row.anno$C))))
  
  # generate gene labels for genes-of-interest
  anno_row <- gtf[rownames(d),'transcript_name']

  hm <- pheatmap::pheatmap(d, scale='none', cluster_cols=F, labels_row = anno_row,  border_color=NA, fontsize_row=4,
                           color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),show_colnames=F,annotation_colors=color.anno,
                           show_rownames=T, treeheight_col=10, treeheight_row=10, clustering_method='ward.D2', annotation_legend=F,
                           angle_col='90', main=paste0(timepoint, 'hr DEG'), 
                           cutree_rows=k, filename=NA, height=4, width=4)
  
  # make gene lists
  res <- read.delim(paste0('Data/03_gm_ora_timepoint_',timepoint,'hrs.txt'))
  res <- res[which((res$overlap > 4) & (res$padj < 0.05) & (res$size <= 600)),]
  res$NonredundantTerms <- T
  res$NonredundantTerms[which(res$pathway_set=='bp')] <- DetectNonRedundantTerms(res$overlapGenes[which(res$pathway_set=='bp')])
  res <- res[which(res$NonredundantTerms),]
  res$overlapGenes <- strsplit(res$overlapGenes, ',')
  p.data <- slice_head(res, by=c('pathway_set','clust'), n=3)
  p.data$label <- paste0(toupper(substr(p.data$pathway,0,1)),  substr(tolower(p.data$pathway),2,80)) %>% stringr::str_wrap(., 30) #, ' (q=', formatC(res$padj, 2, format='e'),')') 
  p.data$pw_ordered <- paste0(rev(formatC(01:nrow(p.data), width=2, flag='0')), '_', p.data$label) %>% factor()
  col <- structure(ggcolor(4), names=c(1:4))
  p1 <- ggdraw() + draw_grob(hm$gtable)
  p3 <- ggplot(p.data[which(p.data$pathway_set=='hallmark'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
    geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='hallmark')]), position='right') + theme_light() +
    scale_size(limits=c(0.01,0.4)) + scale_color_manual(values=col) + 
    labs(title='Pathway', x= 'q value (-log10)', y='', size='PW genes', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
  p4 <- ggplot(p.data[which(p.data$pathway_set=='tf'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
    geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='tf')]), position='right') + theme_light() +
    scale_size(limits=c(0.01,0.4)) + scale_color_manual(values=col) + 
    labs(title='Transcription factor', x= 'q value (-log10)', y='', size='PW genes', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
  return(list(p1, p3, p4))
})
p2 <- plot_grid(plotlist=unlist(p, recursive=F), align='hv', rel_heights=c(1,1,1), rel_widths = c(2,1,1), axis='blr', ncol=3)
ggsave(plot=p2, paste0('Figures/03_bubl_pw_hallmark_v1.4.pdf'), height=15, width=15)


## Make new pw dotplots
timepoints <- c(0,6,24,72)
p <- lapply(timepoints, function(timepoint){
  res <- read.delim(paste0('Data/03_gm_ora_timepoint_',timepoint,'hrs.txt'))
  
  # make gene lists
  res <- res[which((res$overlap > 4) & (res$padj < 0.05) & (res$size <= 600)),]
  res$NonredundantTerms <- T
  res$NonredundantTerms[which(res$pathway_set=='bp')] <- DetectNonRedundantTerms(res$overlapGenes[which(res$pathway_set=='bp')])
  res <- res[which(res$NonredundantTerms),]
  res$overlapGenes <- strsplit(res$overlapGenes, ',')
  p.data <- slice_head(res, by=c('pathway_set','clust'), n=3)
  p.data$label <- paste0(toupper(substr(p.data$pathway,0,1)),  substr(tolower(p.data$pathway),2,80)) %>% stringr::str_wrap(., 30) #, ' (q=', formatC(res$padj, 2, format='e'),')') 
  p.data$pw_ordered <- paste0(rev(formatC(01:nrow(p.data), width=2, flag='0')), '_', p.data$label) %>% factor()
  col <- structure(ggcolor(4), names=c(1:4))
  p3 <- ggplot(p.data[which(p.data$pathway_set=='hallmark'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
    geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='hallmark')]), position='right') + theme_light() +
    scale_size(limits=c(0.01,0.4)) + scale_color_manual(values=col) + 
    labs(title='Pathway', x= 'q value (-log10)', y='', size='PW genes', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
  p4 <- ggplot(p.data[which(p.data$pathway_set=='tf'),], aes(x=-log10(padj), y=pw_ordered, color=factor(clust))) + geom_segment(aes(yend=pw_ordered, xend=0)) + 
    geom_point(aes(size=recall)) + scale_y_discrete(labels= rev(p.data$label[which(p.data$pathway_set=='tf')]), position='right') + theme_light() +
    scale_size(limits=c(0.01,0.4)) + scale_color_manual(values=col) + 
    labs(title='Transcription factor', x= 'q value (-log10)', y='', size='PW genes', color='Cluster') + xlim(c(0, 1.1*max(-log10(p.data$padj))))
  return(list(p3, p4))
})
p2 <- plot_grid(plotlist=unlist(p, recursive=F), align='hv', rel_widths = c(1,1), axis='none', ncol=2)
ggsave(plot=p2, paste0('Figures/03_bubl_pw_hallmark_v1.2.pdf'), height=12, width=6.5)




###########

## Make gene lists
glist <- list(FA_matabolism=pw$hallmark[['HALLMARK_FATTY_ACID_METABOLISM']],
              Glycolysis=pw$hallmark[['HALLMARK_GLYCOLYSIS']],
              OXPHOS=pw$hallmark[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']],
              Apoptosis=pw$hallmark[['HALLMARK_APOPTOSIS']],
              Inflammation=pw$hallmark[['HALLMARK_INFLAMMATORY_RESPONSE']])

## scale and windsorize normalized lcpm by timepoint
timepoints <- c(0,6,24,72)
lcpm <- read.delim('Data/02_lcpm_tx.txt')
d <- lapply(timepoints, function(timepoint){
  cols <- grep(paste0('_',timepoint,'_'),colnames(lcpm))
  y <- t(scale(t(lcpm[, cols])))
  y <- apply(y, 2, function(x) windsorize(x))
  return(y)
})


## Individual lines for each sample
p <- lapply(names(glist), function(pathway_name){
  tx <- gtf$transcript_id[which(gtf$gene_name %in% glist[[pathway_name]])]
  df <- do.call(rbind, lapply(d, function(y) colMeans(y[intersect(tx, rownames(y)),], na.rm=T)))
  colnames(df) <- str_remove(colnames(df), '_0')
  rownames(df) <- timepoints    
  df2 <- reshape2::melt(df)
  df2$group <- factor(substr(df2$Var2, 0, 1), levels=rev(c('G','N','S','D')))
  p.tmp <- ggplot(df2, aes(x=Var1, y=value, group=Var2, fill=group, color=group)) + 
    geom_smooth(se=F, linewidth=0.6) +  scale_x_continuous(breaks=timepoints)+
    theme_light() + labs(x='Time', y='Mean scaled expr.', title=pathway_name)
  return(p.tmp)
})
p <- cowplot::plot_grid(plotlist=p, ncol=2, align='v')
ggsave(plot=p, height=8, width=8, 'Figures/04_pathway_trendlines_allLines.pdf')

## grouped lines for each genotype
p <- lapply(names(glist)[1:4], function(pathway_name){
  tx <- gtf$transcript_id[which(gtf$gene_name %in% glist[[pathway_name]])]
  #tx <- intersect(tx, deg)
  df <- do.call(rbind, lapply(d, function(y) colMeans(y[intersect(tx, rownames(y)),], na.rm=T)))
  #df <- do.call(rbind, lapply(d, function(y) apply(y[intersect(tx, rownames(y)),],2, median))) #function(z) quantile(z, 0.5))))
  colnames(df) <- str_remove(colnames(df), '_0')
  rownames(df) <- timepoints    
  df2 <- reshape2::melt(df)
  df2$group <- factor(substr(df2$Var2, 0, 1), levels=rev(c('G','N','S','D')))
  p.tmp <- ggplot(df2, aes(x=Var1, y=value, fill=group, color=group)) + geom_smooth(alpha=0.4, linewidth=0.7) + 
    geom_smooth(se=F, linewidth=0.7) + scale_x_continuous(breaks=timepoints)+
    theme_light() + labs(x='Time', y='Mean scaled expr.', title=pathway_name)
  return(p.tmp)
})
p <- cowplot::plot_grid(plotlist=p, ncol=2, align='v')
ggsave(plot=p, height=3, width=8, 'Figures/04_pathway_trendlines.pdf')



#####
### Make heatmaps for the trendline pathway plots


## select genes to label
genes2label <- list(FA_matabolism=c('ACAA1-203','ACADVL-201','ACADVL-206','ACADVL-216','ACADVL-236','ACOT2-204','ACAT2-201','ACAT2-202','DHCR24-201','ELOVL5-201','FASN-201'),
                    Glycolysis=c('ALDOA-203','ALDOA-206','ALDOA-214','ALDOA-219','ENO1-201','ENO1-203','ENO1-206','ENO2-201','ENO2-205','G6PD-202','HK2-201','LDHA-201','LDHA-202','LDHA-203','LDHA-206','LDHA-221','LDHA-223','PGAM1-201','PGK1-201','PGK1-203','PKM-202','PKM-212','PKM-214','PKM-216','VEGFA-205','VEGFA-211','VEGFA-212'),
                    OXPHOS=c('ATP5F1B-201','ATP5F1C-202','ATP5F1D-202','ATP5MC2-202','ATP5ME-201','ATP5MC3-203','ATP5PB-206','ATP5PF-201','ATP5PO-201','ATP6AP1-201','ATP6V0B-205','ATP6V0C-201','ATP6V1E1-201','ATP6V1F-201','ATP6V1H-204','BAX-202','BAX-207','BAX-211','COX17-201','COX4I1-201','COX5A-204','COX5B-201','COX6A1-201','COX7C-201','COX7C-201','COX8A-201','CS-201','CS-214','FH-201','IDH1-201','IDH2-201','IDH2-202','IDH2-204','IDH3A-206','IDH3B-201','MDH1-201','MDH1-215','MDH2-201','OGDH-201','OGDH-207','SDHA-213','SDHB-201','SDHD-206'))
deg.label <- do.call(rbind, lapply(1:3, function(i){
    deg <- do.call(rbind, lapply((i+1):4, function(j){
        genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
        genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
        deg <- do.call(rbind, lapply(c(0,6,24,72), function(timepoint){
            de <- read.delim(paste0('Data/02_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))
            return(de[which(de$FDR < 0.01 & abs(de$logFC) > 0),c('FDR','gene_name','transcript_name')])
        }))
        return(deg)
    }))
    return(deg)
})) 
deg.label <- deg.label[order(deg.label$FDR),]
deg.label <- slice_head(deg.label[which(deg.label$transcript_name %in% unlist(genes2label)),], by='gene_name', n=1)[,3] %>% unique()


## Get DE genes to use for the heatmaps
deg <- (lapply(1:3, function(i){
    deg <- (lapply((i+1):4, function(j){
        genotype1 <- c('DKO','NR4A1','SOCS3','GTC')[i]
        genotype2 <- c('DKO','NR4A1','SOCS3','GTC')[j]
        deg <- lapply(c(0,6,24,72), function(timepoint){
            de <- read.delim(paste0('Data/02_de_tx_',timepoint,'hrs_',genotype1,'-',genotype2,'.txt'))
            return(de$transcript_id[which(de$FDR < 0.01 & abs(de$logFC) > 0)])
        }) %>% unlist() %>% unique()
        return(deg)
    }))
    return(deg)
})) 
deg <- unique(unlist(deg))  


## Make heatmaps
#combine the data into one data frame
d <- do.call(cbind, d)
p <- lapply(names(glist)[1:3], function(pathway_name){
    tx <- gtf$transcript_id[which(gtf$gene_name %in% glist[[pathway_name]])]
    tx <- intersect(tx, deg)
    d.tmp <- d[intersect(tx, rownames(d)),]
    d.tmp <- d.tmp[!is.na(rowMeans(d.tmp)),]
    
    col.labels <- gtf[tx, 'transcript_name']
    col.labels[which(!(col.labels %in% intersect(deg.label, genes2label[[pathway_name]])))] <- ''
    sample.order <- as.character(meta$sample[with(meta, order(genotype, time, replicate))])
    
    hm72 <- pheatmap::pheatmap((d.tmp[,48:37]), scale='none', cluster_cols=F, show_rownames=F, #labels_row=col.labels,
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),fontsize_col=5,
                               show_colnames=T, treeheight_col=15, treeheight_row=15, clustering_method='average', annotation_legend=F,
                               angle_col='90', main=pathway_name, annotation_col = meta[colnames(d), 2, drop=F], col_breaks=cumsum(rep(12,3)),
                               cutree_cols=1, filename=NA, height=4, width=4)
    df <- data.frame(y=(gtf[hm72$tree_row$labels[hm72$tree_row$order],4]), labels=(col.labels[(hm72$tree_row$order)]))
    o72=hm72$tree_row$order
    df <- data.frame(y=length(o72):1, label=gtf[rownames(d.tmp)[o72],'transcript_name'] )
    df$label[which(!(df$label %in% intersect(deg.label, genes2label[[pathway_name]])))] <- ''
    hm72.label <- ggplot(df, aes(x=0, y=y, label=label)) + geom_point(shape='-') + 
        ggrepel::geom_text_repel(nudge_x=2, min.segment.length=0, size=1.2, max.overlaps=50, segment.size=0.1) + theme_test()
    hm72 <- ggdraw() + draw_grob(hm72$gtable)
    hm72 <- cowplot::plot_grid(hm72,hm72.label, ncol=2, align='v', rel_widths=c(1,0.4)) 
    
    hm00 <- pheatmap::pheatmap((d.tmp[,1:12]), scale='none', cluster_cols=F,  show_rownames=F,
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100),fontsize_col=5,
                               show_colnames=T, treeheight_col=15, treeheight_row=15, clustering_method='average', annotation_legend=F,
                               angle_col='90', main=pathway_name, annotation_col = meta[colnames(d), 2, drop=F], col_breaks=cumsum(rep(12,3)),
                               cutree_cols=1, filename=NA, height=4, width=4)
    df <- data.frame(y=(gtf[hm00$tree_row$labels[hm00$tree_row$order],4]), labels=(col.labels[(hm00$tree_row$order)]))
    o00=hm00$tree_row$order
    df <- data.frame(y=length(o00):1, label=gtf[rownames(d.tmp)[o00],'transcript_name'] )
    df$label[which(!(df$label %in% intersect(deg.label, genes2label[[pathway_name]])))] <- ''
    hm00.label <- ggplot(df, aes(x=0, y=y, label=label)) + geom_point(shape='-') + 
        ggrepel::geom_text_repel(nudge_x=2, min.segment.length=0, size=1.2, max.overlaps=50, segment.size=0.1) + theme_test()
    hm00 <- ggdraw() + draw_grob(hm00$gtable)
    hm00 <- cowplot::plot_grid(hm00,hm00.label, ncol=2, align='v', rel_widths=c(1,0.4)) 
    
    return(list(hm00, hm72))
})
p <- unlist(p, F)
ggsave(plot=cowplot::plot_grid(plotlist=p, ncol=6, align='v', axis='lbt'), 
       height=3, width=24, 'Figures/04_pathway_heatmaps_v1.4.6.pdf')







  
  
