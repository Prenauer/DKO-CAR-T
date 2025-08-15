

combineMatrices <- function(m1, m2, suffix='-1'){
    g1 <- setdiff(rownames(m1), rownames(m2))
    g2 <- setdiff(rownames(m2), rownames(m1))
    
    m1 <- rbind(m1, matrix(0, nrow=length(g2),ncol=ncol(m1), 
                           dimnames=list(g2,colnames(m1))))
    m2 <- rbind(m2, matrix(0, nrow=length(g1),ncol=ncol(m2), 
                           dimnames=list(g1,colnames(m2))))
    m2 <- m2[rownames(m1),]
    colnames(m2) <- stringr::str_replace(colnames(m2), '-1', suffix)
    m1 <- cbind(m1,m2)
    return(m1)
}
DetermineOptimalClusters <- function(so, resolution=NULL, plot.only=F, graph.name='SCT_snn',
                                     plot.res=NULL, optimalResolution=NULL, seed=42){
    if(!plot.only){
        ## Optimize clustering (calculating WSS and silhouette-width)
        plan('sequential')
        xy <- so@reductions$umap@cell.embeddings %>% data.frame()
        d <- dist(xy)#Matrix::Matrix(so[[graph.name]])
        CalcWSS <- function(x,y) (length(x)-1)*(var(x)+var(y))
        optimalResolution <- do.call(rbind, lapply(resolution, function(res){
            tmp <- FindClusters(so, graph.name = graph.name, algorithm = 4, 
                                resolution = res, random.seed = seed)
            wss <- reframe(cbind(xy,cluster=Idents(tmp)), .by='cluster', 
                           wss=CalcWSS(umap_1,umap_2))$wss %>% sum()
            asw <- mean(cluster::silhouette(as.integer(Idents(tmp)), dist=d)[,3])
            cat(res,'\n')
            return(data.frame(res=res, wss=wss, asw=asw, n_clust=length(levels(Idents(tmp)))))
        }))
    }
    
    ## plot cluster-n, WSS, and ASW as a function of resolution
    if(is.null(plot.res)){
        res <- optimalResolution[with(optimalResolution, 
                                      order(n_clust,-asw,wss)),]
        r <- reframe(res, .by='n_clust', 
                     asw.pick=head(res[which(asw==max(asw))],1))[,2] %>% unlist()
    }
    else{
        r <- plot.res
    }
    
    p.resCheck <- cowplot::plot_grid(ggplot(optimalResolution, aes(x=res, y=n_clust)) + 
                                         geom_point(color='tomato2') + 
                                         geom_path(color='tomato2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('# of clusters',20)) + 
                                         theme_test(),
                                     ggplot(optimalResolution, aes(x=res, y=wss)) + 
                                         geom_point(color='cyan2') + 
                                         geom_path(color='cyan2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('Within-cluster sum-of-squares (wss)',20)) +
                                         theme_test(),
                                     ggplot(optimalResolution, aes(x=res, y=asw)) + 
                                         geom_point(color='orange2') + 
                                         geom_path(color='orange2') + 
                                         geom_vline(xintercept=r, linetype='dashed', 
                                                    color='gray25', alpha=0.75) +
                                         labs(x='Clustering resolution', 
                                              y=stringr::str_wrap('Average silhouette width',20)) + 
                                         theme_test(),
                                     ncol=1, align='hv')
    cat('Optimal res = ', paste0(sort(r),collapse=','),'\n')
    
    return(list(res=optimalResolution, plot=p.resCheck))
}


## Make function for generating a customized Dot Plot
customDotPlot <- function(so, genes, group.order=NULL, 
                          group.by='ct', scale.max=60){
    ## Reorder the cells according to the provided order
    if(!is.null(group.order)){
        so@meta.data[,group.by] <- factor(so@meta.data[,group.by], 
                                          levels=group.order)
    }
    ## Make Dotplot of the provided genes
    DotPlot(so,genes, 
            group.by=group.by, scale.max=scale.max) + 
        # remove labels
        labs(x=NULL, y=NULL) +
        # customize legend presentation
        theme(legend.position='bottom', legend.key.spacing=unit(0.5,'lines'), 
              legend.key.height=unit(0.5,'lines'), 
              legend.key.width=unit(0.5,'lines'), 
              legend.text=element_text(size=7), 
              legend.key.spacing.y=unit(0.01,'lines'),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              legend.title=element_text(size=10)) +
        # change dot color-scheme
        scale_color_viridis_c(option='inferno', direction=-1, end=0.9) +
        # set order of which lengend is shown first
        guides(size = guide_legend(nrow = 2), alpha = guide_legend(nrow = 2))
}

plotVolcano <- function(de, lfc.thresh=1, p.thresh=0.01, pt.size=2, title=NULL,
                        text.size=2.2, rastr=T, xlim.expansion=1.3){
    de$p <- de[,'padj']
    de$lfc <- de[,'logFC']
    de$feature <- de[,'feature']
    de <- de[order(de$pval),]
    genes2label <- c(head(de$feature[which((de$lfc > lfc.thresh)&(de$p < p.thresh))],5),
                     head(de$feature[which((de$lfc < -lfc.thresh)&(de$p < p.thresh))],5))
    de <- de[order(abs(de$lfc), decreasing=T),]
    genes2label <- c(genes2label,
                     head(de$feature[which((de$lfc > lfc.thresh)&(de$p < p.thresh))],5),
                     head(de$feature[which((de$lfc < -lfc.thresh)&(de$p < p.thresh))],5))
    genes2label <- unique(genes2label)
    
    # set min pval limit for plotting
    plot.min.p <- min(de$p[which(de$p != 0)], na.rm=T)*0.1
    de$p[which(de$p == 0)] <- plot.min.p
    
    # get plot size limits
    plot.xlim <- xlim.expansion*range(de$lfc, na.rm=T)
    plot.ylim <- max(-log10(de$p), na.rm=T)*c(-0.02,1.2)
    
    p <- ggplot(de, aes(x=lfc, y=-log10(p))) 
    if(rastr){
        p <- p + geom_point_rast(color='gray20', alpha=1, stroke=0, size=pt.size*1.5) +
            geom_point_rast(color='gray', alpha=1, stroke=0, size=pt.size) +
            geom_point_rast(data=de[which(de$lfc > lfc.thresh & de$p < p.thresh),],
                            color='firebrick', alpha=1, stroke=0, size=pt.size) +
            geom_point_rast(data=de[which(de$lfc < -lfc.thresh & de$p < p.thresh),],
                            color='steelblue', alpha=1, stroke=0, size=pt.size)
    }
    if(!rastr){
        p <- p + geom_point(color='gray', alpha=0.5, stroke=0, size=pt.size) +
            geom_point(data=de[which(de$lfc > lfc.thresh & de$p < p.thresh),],
                       color='firebrick', alpha=0.8, stroke=0, size=pt.size) +
            geom_point(data=de[which(de$lfc < -lfc.thresh & de$p < p.thresh),],
                       color='steelblue', alpha=0.8, stroke=0, size=pt.size)
    }
    p <- p + ggrepel::geom_text_repel(data=de[which(de$lfc < 0 & de$feature %in% genes2label),],
                                      aes(label=feature), min.segment.length=0,
                                      nudge_x=-0.5, nudge_y=0.5, max.iter = 1e5,
                                      segment.alpha=0.7, segment.size=0.01,
                                      size=text.size, max.overlaps=50, fontface='italic') +
        ggrepel::geom_text_repel(data=de[which(de$lfc > 0 & de$feature %in% genes2label),],
                                 aes(label=feature), min.segment.length=0,
                                 nudge_x=0.5, nudge_y=0.5, max.iter = 1e5,
                                 segment.alpha=0.7, segment.size=0.01,
                                 size=text.size, max.overlaps=50, fontface='italic') +
        xlim(plot.xlim) + ylim(plot.ylim) +
        theme_classic() + labs(x='Log2 fold-change', y='Adj. p (-log10)', title=title)
    return(p)
}

RunNetPW <- function(pw, plot.title=NULL, n.pw.per.group=5, seed=42){
    # Select top pathways from each comparison
    top.pw <- slice_head(pw, by='group', n=n.pw.per.group)$pathway
    # Subset by top pathways
    pw <- pw[which(pw$pathway %in% top.pw),]
    # Generate edges 
    e <- rbind(
        # from celltype-group to comparison
        data.frame(from=pw$ct, to=pw$group, size=0.1),
        # from comparison to pathway
        data.frame(from=pw$group, to=pw$pathway, size=pw$NES/max(pw$NES, T)))
    # Generate vertices
    v <- rbind(data.frame(name=pw$pathway,type='pathway'),
               data.frame(name=pw$ct,type='celltype'),
               data.frame(name=pw$group,type='comparison'))
    v <- unique(v)
    
    
    # Create igraph object
    ig <- graph_from_data_frame( e, vertices=v )
    # generate layout
    set.seed(seed)
    ig2 <- (ig %>% add_layout_(with_fr()))
    # convert to ggnetwork
    net <- ggnetwork(ig2)
    # create column that states the type of node
    net$class <- 'edge'
    net$class[which(net$x==net$xend & net$y==net$yend)] <- 'node'
    
    # Reformat pathway names
    net$name[which(net$type=='pathway')] <- net$name[which(net$type=='pathway')] %>%
        str_remove(.,'GOBP_') %>% str_remove(.,'REACTOME_') %>%
        str_remove(.,'KEGG_') %>% str_remove(.,'WP_') %>% 
        str_remove(.,'BIOCARTA_') %>% str_replace_all(., '_', ' ')
    
    # exclude celltype nodes and edges from the visualization; 
    #   only show those for the comparisons.
    net <- net[which(net$type != 'celltype'),]
    
    # adjust names of 'comparison' nodes
    net$node_label=NA
    net$node_label[which(net$class=='node')] <- 
        str_split_i(net$name[which(net$class=='node')], ': ',1) 
    
    ## Make plot
    # prep colors
    cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                      names=c('GTC','SOCS3','NR4A1','DKO'))
    legend.labs <- structure(c('GTC','SOCS3','NR4A1','DKO'),
                             names=c('GTC','SOCS3','NR4A1','DKO'))
    # fix names
    net$gp2 <- str_split_i(net$name, '-', 2)   
    # make plot
    ggplot(net, aes(x = x, y = y)) + 
        geom_edges(color="#DDDDDD",
                   aes(linewidth=size, xend = xend, yend = yend)) +
        geom_edges(data=net[which(net$class=='edge' & net$type=='comparison'),], 
                   arrow = arrow(length = unit(0.25, 'lines'), 
                                 type = "closed"), #curvature=0.4, 
                   aes(linewidth=size, color=gp2, xend = xend, yend = yend)) +
        scale_linewidth(range=c(0.01,1)) +
        theme_blank() +
        geom_nodes(data=net[which(net$class=='node' & net$type=='comparison'),], 
                   shape=22, color='gray30', size=6,
                   aes(fill=gp2)) + 
        geom_nodes(data=net[which(net$class=='node' & net$type=='pathway'),], 
                   shape=21, color='gray30',fill='#20854EFF', size=3) + 
        geom_nodetext_repel(data=net[which(net$class=='node' & net$type=='comparison'),], 
                            aes(label=node_label), size=3, min.segment.length=0, 
                            nudge_y=0.5, fontface='bold',
                            segment.size=0.1, segment.alpha=0.5, max.overlaps=50) + 
        geom_nodetext_repel(data=net[which(net$class=='node' & net$type=='pathway'),], 
                            aes(label=str_wrap(name,20)), size=1.7, min.segment.length=0, 
                            nudge_y=0.5, lineheight=0.8,
                            segment.size=0.1, segment.alpha=0.5, max.overlaps=50) + 
        scale_fill_manual(values=cols, labels=legend.labs) + 
        scale_color_manual(values=cols, labels=legend.labs) + 
        labs(title=plot.title, fill='DE', 
             linewidth=str_wrap('Pathway enrichment score',8)) + 
        guides(fill=guide_legend(order = 1),linewidth=guide_legend(order=2), 
               color=F) +
        theme(legend.key.height=unit(0.2,'lines'), legend.spacing.x=unit(0.2,'lines')) 
}


MakeDotPlot <- function(g, title='', thresh=1/5, method='ward.D2'){
    require(ggdendroplot)
    g <- intersect(hvg, g)
    m <- t(so$RNA$data[g,])
    m <- splitAsList(m, factor(so$id, levels=unique(so$id)))
    m <- do.call(rbind, lapply(m, colMeans))
    rownames(m) <- str_replace(rownames(m), '_','-')
    res.pca <- prcomp(m, scale=T)
    
    g <- g[(apply(so$RNA$data[g,], 1, PercentAbove, threshold=0) > thresh)]
    
    m <- t(so$RNA$data[g,])
    m <- splitAsList(m, factor(so$id, levels=unique(so$id)))
    #m <- do.call(rbind, lapply(m, function(x){apply(x, 2, PercentAbove, threshold=0)}))
    m <- do.call(rbind, lapply(m, colMeans)) %>% t()
    #rownames(m) <- unique(so$id)
    cl.x <- hclust(dist(m), method=method)
    cl.y <- hclust(dist(t(m)), method=method)
    so$id <- factor(so$id, levels=cl.y$labels[cl.y$order])
    g <- g[cl.x$order]
    
    p1 <- DotPlot(so,g,group.by='id') + #coord_flip() + 
        labs(x=NULL, y=NULL, title=title) +
        theme(legend.position='bottom', legend.key.spacing=unit(0.5,'lines'), 
              legend.key.height=unit(0.5,'lines'), legend.key.width=unit(0.5,'lines'), 
              legend.text=element_text(size=7), legend.key.spacing.y=unit(0.01,'lines'),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'),
              legend.title=element_text(size=10)) +
        scale_color_distiller(palette = "RdYlBu") +
        guides(size = guide_legend(nrow = 2), alpha = guide_legend(nrow = 2))
    p1.dend <- ggplot() + geom_dendro(cl.y, pointing='side') + theme_void()
    d <- data.frame(res.pca$x[,1:2], id=rownames(res.pca$x),
                    ct=str_split_i(rownames(res.pca$x),'-',1),
                    geno=str_split_i(rownames(res.pca$x),'-',2))
    d$ct <- str_replace(d$ct, 'Tcell', 'Act.T')
    p2 <- ggplot(d, aes(x=PC1, y=PC2, color=geno)) + geom_point(size=4) + 
        theme_bw() + theme(legend.position='none') + 
        scale_color_manual(values=
                               structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                                         names=c('GTC','SOCS3','NR4A1','DKO'))) +
        ggrepel::geom_text_repel(aes(label=ct), size=4, min.segment.length=0) + 
        labs(title=title)
    #p <- cowplot::plot_grid(p1,p2, align='h', axis='lbt',ncol=2, rel_widths=c(0.6,0.4))
    p <- list(p1,p1.dend, p2)
    return(p)
}


GamPlot <- function(g,pheno=NULL,title=NULL,celltypes= NULL, raster=T){ 
    if(is.null(celltypes)) celltypes <- unique(so$groups)
    if(g %in% rownames(so$es.tf)){
        d <-data.frame(y=e[,pheno], x=so$es.tf$data[g,], pheno=pheno,
                       ct=so$groups, gp=so$sample,ct.sub=so$ct)
    }
    if(g %in% rownames(so$RNA)){
        d <-data.frame(y=e[,pheno], x=so$RNA$data[g,], pheno=pheno,ct=so$groups,
                       ct.sub=so$ct, umi=so$nCount_RNA, gp=so$sample)
    }
    d <- d[which(d$ct %in% celltypes),]
    d <- d[which(d$x > 0 & d$y > 0),]
    cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                      names=c('GTC','SOCS3','NR4A1','DKO'))
    d$gp <- factor(d$gp, levels=c('GTC','SOCS3','NR4A1','DKO'))
    smooth.x <- seq(quantile(d$x, 0.01),quantile(d$x, 0.99), 
                    diff(quantile(d$x, c(0.01,0.99)))/100)
    
    if(is.null(title)) title <- paste0(pheno,' signature ~ ',g)
    p <- ggplot(d[order(d$gp),], aes(y=y,x=x)) 
    if(raster) p <- p + ggrastr::geom_point_rast(size=0.01, aes(color=gp))
    if(!raster) p <- p + geom_point(size=0.01, aes(color=gp))
    p <- p +
        theme_classic() + 
        scale_color_manual(values=cols) +
        geom_smooth(method.args=list(family='Gamma', method='ML', select=T),
                    method='gam', formula=y~ x + s(x, bs="bs"), alpha=0.1, 
                    fullrange=F,xseq=smooth.x) +
        labs(x=paste0(g,' signature'), y=paste0(str_replace(pheno,'-TARGET-GENES','_TF'),' signature'),
             title=title)+
        theme(legend.position='none', title=element_text(size=7))
    
    p <- ggExtra::ggMarginal(p, groupColour=T, margins='both', groupFill=T, size=4, alpha=0.2) 
    return(p)
}


GamEffectPlotTF <- function(x=NULL, y, resid=NULL, title=NULL, raster=T){ # note make sure "input" is available
    ## if input is a TF signature
    if(x %in% rownames(so$es.tf)){
        d <-data.frame(y=e[,y], x=so$es.tf$data[x,], 
                       ct= factor(so$ct),
                       s=factor(so$sample, levels=c('GTC','NR4A1','SOCS3','DKO')))
        x <- str_replace(x, '-TARGET-GENES', '_TF')
        lab.y <- paste0(y,' signature')
        lab.x <- paste0(x, ' signature')
        title <- paste0(y, '~', x, ' relationship')
        d <- d[which(d$x > 0 & d$y > 0),]
        fit <- mgcv::bam(y ~ ct + s + s(x, bs="ps",k=10),  
                         control=list(keepData=F),
                         data=d,family=Gamma(link='log'), method='fREML')
    }
    
    ## Get residuals
    pr <-  cbind(predict(fit, newdata=d, type='terms'), resid=residuals(fit))
    if(is.null(resid)) resid <- 'x'
    if(resid=='x'){
        d$p <- rowSums(pr[,c('s(x)','resid')])[rownames(d)]
    }
    if(resid=='both'){
        d$p <- rowSums(pr[,c('s(x)','s','resid')])[rownames(d)]
    }
    
    ## create annotation for stats
    anno <- structure(c(cor.test(d$x,d$p)[c('estimate','p.value')], summary(fit)$r.sq),
                      names=c('cor_rsq','pval','fit_rsq'))
    #if(sign(cor.test(d$x,d$p)$estimate
    anno <- paste0('Correlation:\nr2 = ',formatC(anno[[1]], digits=3),
                   '\nadj.p = ',formatC(anno[[2]], digits =2,format='e'))
    cat(anno)
    
    smooth.x <- seq(quantile(d$x, 0.01),quantile(d$x, 0.99), 
                    diff(quantile(d$x, c(0.01,0.99)))/100)
    anno.x <- min(d$x) + (diff(range(d$x))*0.8)
    anno.y <- min(d$p) + (diff(range(d$p))*0.25)
    
    ## Make plot
    cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                      names=c('GTC','SOCS3','NR4A1','DKO'))
    #p <- ggplot(d[order(d$s),], aes(y=p,x=x))
    set.seed(1)
    p <- ggplot(d[sample(1:nrow(d),nrow(d),F),], aes(y=p,x=x))
    if(!raster) p <- p + geom_point(aes(color=s)) 
    if(raster) p <- p + ggrastr::geom_point_rast(size=0.01, aes(color=s))
    p <- p + theme_classic() + 
        scale_color_manual(values=cols) +
        #scale_fill_brewer(palette='reds') +
        geom_smooth(method='gam', formula=y~ s(x, bs="ps",k=10), alpha=0.1, 
                    color='gray20', fullrange=F,xseq=smooth.x) +
        labs(x=lab.x, y=lab.y, title=title)+
        annotate('text',label=anno, x=anno.x, y=anno.y, size=2.5,
                 lineheight=0.7) +
        theme(legend.position='none', title=element_text(size=7))
    
    p <- ggExtra::ggMarginal(p, groupColour=T, margins='both', groupFill=T, size=4, alpha=0.2) 
    return(p)
    
}
GamPlotTF <- function(x,y=NULL,title=NULL,celltypes= NULL){ # note make sure "input" is available
    if(is.null(celltypes)) celltypes <- unique(so$groups)
    d <-data.frame(y=e[,y], x=so$es.tf$data[x,], pheno=y,ct=so$groups,
                   ct.sub=so$ct,
                   #sf=log(so$nCount_RNA),
                   gp=so$sample)
    if(sum(celltypes %in% d$ct.sub) > sum(celltypes %in% d$ct)) d$ct <- d$ct.sub
    d <- d[which(d$ct %in% celltypes),]
    d <- d[which(d$x > 0 & d$y > 0 ),]
    cols <- structure(c('#B2B2B2','#E9AC4C','#496AB4','#D85446'),
                      names=c('GTC','SOCS3','NR4A1','DKO'))
    d$gp <- factor(d$gp, levels=c('GTC','SOCS3','NR4A1','DKO'))
    smooth.x <- seq(quantile(d$x, 0.01),quantile(d$x, 0.99), 
                    diff(quantile(d$x, c(0.01,0.99)))/100)
    
    if(is.null(title)) title <- paste0(y,' ~ ',
                                       str_replace(x,'-TARGET-GENES','_TF'), 
                                       ' relationship')
    p <- ggplot(d[order(d$gp),], aes(y=y,x=x)) +
        
        ggrastr::geom_point_rast(size=0.01, aes(color=gp))+
        #geom_density_2d(data=d[which(d$gp=='DKO'),], color='#D85446', bins=6) +
        theme_classic() + 
        scale_color_manual(values=cols) +
        #scale_fill_brewer(palette='reds') +
        geom_smooth(method.args=list(family=Gamma(link='log'), method='REML', select=T),
                    method='gam', formula=y~ s(x, bs="ps",k=10), alpha=0.1, 
                    fullrange=F,xseq=smooth.x, color='gray20') +
        labs(y=paste0(y,' signature'), x=paste0(str_replace(x,'-TARGET-GENES','_TF'),' signature'),
             title=title)+
        theme(legend.position='none', title=element_text(size=7))
    
    p <- ggExtra::ggMarginal(p, groupColour=T, margins='both', groupFill=T, size=4, alpha=0.2) 
    return(p)
}


