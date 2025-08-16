library(BiocManager)
library(edgeR)
library(stringr)
library(dplyr)
require(rtracklayer)
require(rhdf5)
require(RColorBrewer)
require(ggrepel)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(reshape2)
library(factoextra)
library(cluster)
library(decoupleR)
source('Code/samba_v2.0.R')
ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]


## load data
sgres <- read.delim('Data/dko_invivoScr_GuideLevelResults.txt')
gres <- read.delim('Data/dko_invivoScr_GeneLevelResults.txt')
gres$z_pos <- scale(gres$score_pos)[,1]
write.table(gres, 'Data/dko_invivoScr_GeneLevelResults_v2.txt', sep='\t', quote=F, row.names=F)


## setup dataframe for epistasis analyses
epi_input <- data.frame(Gene = gres$Gene, 
                        GeneA = stringr::str_split(gres$Gene, '_', simplify = T)[,1], 
                        GeneB = stringr::str_split(gres$Gene, '_', simplify = T)[,2], 
                        Samba_qval = gres$qval_pos,
                        dkoAB = gres$z_pos)
sko <- epi_input[which(epi_input$GeneB == '') ,c('GeneA','dkoAB')]
epi_input$skoA = sko$dkoAB[match(epi_input$GeneA, sko$GeneA)]
epi_input$skoB = sko$dkoAB[match(epi_input$GeneB, sko$GeneA)]
epi_input[is.na(epi_input)] <- 0


## DKO score grid
rev.input <- epi_input[,c(1,3,2,4,5,6,7)]
colnames(rev.input) <- colnames(epi_input)#[1:5]
d <- rbind(epi_input, rev.input)
d$epi_score <- (d$dkoAB) - ((d$skoA + d$skoB))
for(i in c(2:3)) d[which(d[,i]==''), i] <- 'GTC'

# heatmap for synergy score
d2 <- reshape2::dcast(unique(d), GeneA ~ GeneB, value.var='epi_score', fill=0)
d2 <- data.frame(d2[,-1], row.names=d2[,1])
windorize=function(x, thresh=2) {
    x[which(x > thresh)] <- thresh
    x[which(x < -thresh)] <- -thresh
    return(x)
}
d2 <- apply(d2, 2, windorize, thresh=2)
pheatmap::pheatmap(d2, scale='none', clustering_method='ward.D2', cutree_row=4, cutree_col=4, angle_col='90',
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")))(100),
                   main='Screen synergy across gene pairs',
                   treeheight_col=20, treeheight_row=20, height=10, width=11, filename='Figs/invivo_scr_chessmap_SynergyScore_v3.7.pdf')

# heatmap for samba score
d2 <- reshape2::dcast(unique(d), GeneA ~ GeneB, value.var='dkoAB', fill=0)
d2 <- data.frame(d2[,-1], row.names=d2[,1])
d2 <- apply(d2, 2, windorize, thresh=1)
pheatmap::pheatmap(d2, scale='none', clustering_method='ward.D2', cutree_row=4, cutree_col=4, angle_col='90',
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")))(100),
                   main='Screen enrichment across gene pairs',
                   treeheight_col=20, treeheight_row=20, height=10, width=11, filename='Figs/invivo_scr_chessmap_sambaScore_v3.7.pdf')


## run epistasis
epi <- data.frame(Gene = epi_input$Gene,
                  Epi_score = (epi_input$dkoAB) - ((epi_input$skoA + epi_input$skoB)),
                  Add_score = (epi_input$skoA + epi_input$skoB))
epi <- epi[grep('_',epi$Gene),]
epi$Epi_pval <- 2*pnorm(abs(epi$Epi_score), lower.tail = F)
epi$Epi_qval <- p.adjust(epi$Epi_pval)


## write results to file
df = merge(epi_input, epi, by = 'Gene')
colnames(df) <- c('Gene','Gene1','Gene2','Samba_qval','DKO_score','Gene1_score','Gene2_score','Epistasis_score','Additive_score','Epistasis_pval','Epistasis_qval')
write.table(df, 'Data/dkot_invivo_epistasis_v3.7.txt', sep = '\t', row.names = F, quote = F)


## plot epistasis results
df <- read.delim('Data/dkot_invivo_epistasis_v3.7.txt')
samba.low <- -2.84
samba.high <- 2.84
epi.low <- -4.26
epi.high <- 4.28
df <- df[order(df$Epistasis_score, decreasing=T),]
tophits = c(df$Gene[(df$DKO_score <  samba.high) & (df$Epistasis_score >  epi.high)] %>% head(4),
            df$Gene[(df$DKO_score >  samba.high) & (df$Epistasis_score >  epi.high)] %>% head(4),
            df$Gene[(df$DKO_score > samba.low) & (df$Epistasis_score < epi.low)] %>% tail(4),
            df$Gene[(df$DKO_score < samba.low) & (df$Epistasis_score < epi.low)] %>% tail(4))
df <- df[order(df$DKO_score, decreasing=T),]
tophits = c(tophits, head(gres$Gene,8),
            df$Gene[(df$DKO_score >  samba.high) & (df$Epistasis_score >  epi.high)] %>% head(4), 
            df$Gene[(df$DKO_score >  samba.high) & (df$Epistasis_score <  epi.high)] %>% head(4),
            df$Gene[(df$DKO_score < samba.low) & (df$Epistasis_score > epi.low)] %>% tail(4),
            df$Gene[(df$DKO_score < samba.low) & (df$Epistasis_score < epi.low)] %>% tail(4))
tophits <- unique(tophits)
df.label <- df[which(df$Gene %in% tophits),]
df.label$textsize <- 1
df.label$textsize[which(df.label$Gene %in% df[(df$Samba_qval < 0.05) & (df$Epistasis_qval < 0.05),'Gene'])] <- 2
df.label.up <- df.label[which(df.label$Epistasis_score > 0),]
df.label.dn <- df.label[which(df.label$Epistasis_score < 0),]
anno.equation <- data.frame(DKO_score= -4.9, Epistasis_score = 4.3, label='Epis = DKO - (SKO1 + SKO2)')
p <- ggplot(df, aes(x=DKO_score, y=Epistasis_score)) + xlim(c(-5,5)) + #ylim(c(-4.5,4)) + 
    rasterize(ggpointdensity::geom_pointdensity(show.legend = F)) + 
    scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
    geom_point_rast(data = df[(abs(df$DKO_score) > samba.high) | (df$Epistasis_qval < 0.05),], size = 2.5, color = 'gray10') + 
    geom_point_rast(data = df[(abs(df$DKO_score) > samba.high) | (df$Epistasis_qval < 0.05),], size = 2, color = 'white') + 
    geom_point_rast(data = df[which(abs(df$DKO_score) > samba.high),], size = 2, color = 'steelblue', alpha=1)+ 
    geom_point_rast(data = df[which(df$Epistasis_qval < 0.05),], size = 2, color = 'firebrick', alpha=1) + 
    geom_point_rast(data = df[(abs(df$DKO_score) > samba.high) & (df$Epistasis_qval < 0.05),], size = 2, color = 'purple', alpha=1) + 
    geom_hline(yintercept = c(epi.low,epi.high), color = 'gray50', linetype = 'dashed') + 
    geom_vline(xintercept = c(samba.low, samba.high), color = 'gray50', linetype = 'dashed') + 
    geom_text(data = anno.equation, aes(label=label), size=3, hjust=0) + 
    labs(x = 'Samba gene-pair z score', y = 'Epistasis z score', title = 'Epistasis') +
    ggnewscale::new_scale('size')  +
    geom_text_repel(data = df.label.up, aes(label = Gene, size=(textsize)), segment.alpha = 0.5,
                    max.iter = 1000000, nudge_x=0.5, nudge_y=0.5,
                    min.segment.length = 0, segment.size = 0.15, max.overlaps = 10) + 
    geom_text_repel(data = df.label.dn, aes(label = Gene, size=(textsize)), segment.alpha = 0.5,
                    max.iter = 1000000, nudge_x=-0.5, nudge_y=-0.5, #fontface='italic', 
                    min.segment.length = 0, segment.size = 0.15, max.overlaps = 10) + 
    scale_size(range=c(2,2.6), guide='none') + 
    cowplot::theme_cowplot() + theme(legend.position = 'none')
ggsave(plot = p, filename = 'Figs/invivo_xy_Epistasis_v3.7.pdf', width = 4, height = 4)







