#' Samba
#'
#' This is an all-in-one function to preprocesses and analyze read count data from a CRISPR
#' screen. The structure of the input data is assumed to have columns labeled "sgRNA" and
#' "Gene", followed by read counts from each sample in appropriately labeled columns.
#'
#' @param data data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param ntc.as.null.dist logical, indicating whether to ese NTCs or all guides to generate a null distribution for gene-level analyses.
#' @param score.method character, indicating whether to use the 'MetaAnalysis' or 'GeneScore' method to perform gene-level analysis.
#' @param test.method character, indicating whether to use 'QLF' or 'LRT' tests.
#' @param GuideMap data.frame that maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.
#' @param file.prefix Prefix for writing output files.
#' @return DGEList object of the edgeR package.
#' @export
Samba <- function(data, design, coefficient = NULL, contrast = NULL, ntc.as.null.dist = T, normalization.method = 'TMMwsp',
                  score.method = 'MetaAnalysis', test.method = 'QLF', GuideMap = NULL, file.prefix = NULL){

    dge <- Preprocess_Samba(data = data, design = design, normalization.method = normalization.method)
    sgRes <- Analyze_Samba_Guides(dge = dge, coefficient = coefficient, method = test.method, contrast = contrast, file.prefix = file.prefix)
    geneRes <- Analyze_Samba_Genes(sgRes = sgRes, ntc.as.null.dist = ntc.as.null.dist, score.method = score.method, file.prefix = file.prefix)

    return(list(GuideResults = sgRes, GeneResults = geneRes))
}


#' Preprocess_Samba
#'
#' This function preprocesses read count data from CRISPR screen data. The structure of the
#' input data is assumed to have columns labeled "sgRNA" and "Gene", followed by read counts
#' from each sample in appropriately labeled columns.
#'
#' @param data data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @return DGEList object of the edgeR package.
#' @export
Preprocess_Samba <- function(data, design, min.guides = 1, pseudocount = 4, group = NULL, normalization.method = 'TMMwsp'){
    cat('Starting preprocessing.\n')
    # Default settings
    hdr.gene = 'Gene'
    hdr.guide = 'sgRNA'
    dispersion.method = 'GLM'

    # prep input
    sgInfo <- data.frame(data[,1:2])
    data <- data.frame(data[,-c(1:2)])
    rownames(data) <- sgInfo[, hdr.guide]
    rownames(sgInfo) <- sgInfo[, hdr.guide]
    design <- as.matrix(design)
  
    # check input
    if (sum(colnames(sgInfo) %in% c(hdr.gene, hdr.guide)) != 2) 
        warning("Count data is missing the columns: \"Guide\" and \"Gene\"!")
    if (nrow(design) != ncol(data)) 
        warning("Design matrix has incorrect number of rows!")
    if ((mode(design) != "numeric") & (mode(design) != "integer")) 
        warning("Design matrix does not have integer/numeric format!")

    # get sample groups
    if(is.null(group)){
        if(ncol(design) == 2) {
            group <- as.integer(design[, 2])
        } else {
            group <- as.integer(design[, ncol(design)])
        }
    }
    if(length(unique(group)) >= (nrow(design))) group = NULL


    # filter data
    if (min.guides > 0) {
	    if (!is.null(group)) {
	        samples.screen <- which(group != 0)
	    } else {
	        samples.screen <- colnames(data)
	    }
	    keep.exprs = rowSums(data[,samples.screen])
	    keep.exprs <- keep.exprs[keep.exprs > min.guides] %>% names()
	    cat(paste0("   Removed ", nrow(data) - length(keep.exprs), 
	        " of ", nrow(data), " guides.\n"))
    } else {
    	keep.exprs <- rownames(data)
    }
    data <- data[keep.exprs,]
    sgInfo <- sgInfo[keep.exprs,]

    # normalize data
    cat("Performing normalization: ",normalization.method,"\n")
    if(normalization.method == 'uq2'){
      data <- UpperQuantNorm2Step(data.frame(data), 100)
      normalization.method <- 'none'
    }
    if(normalization.method == 'total'){
        data <- apply(data, 2, function(x) 1e6*x/sum(x))
        normalization.method <- 'none'
    }
    dge <- edgeR::DGEList(data, remove.zeros = T, group = group)
    dge <- edgeR::calcNormFactors(dge, method = normalization.method)
 
    # generate guide weights based on detection across samples
    #dge <- WeighGuides(dge = dge)

    # calculate size factors and dispersion
    cat('Calculating dispersions ...\n')
    dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
    cat('   Calculated common dispersion\n')
    dge <- edgeR::estimateGLMTrendedDisp(dge, design = design, method = 'bin.loess', span = 1/3)
    cat('   Calculated trended dispersion\n')
    dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
    cat('   Calculated sgRNA-wise dispersion\n')

    # return DGE object
    dge$Design <- design
    dge$GuideMap <- sgInfo
    cat('Preprocessing complete!\n')
    return(dge)
}

#' Analyze_Samba_Guides
#'
#' This function performs guide-level screen enrichment analysis of SAMBA, using a pre-processed DGEList object as an input.
#'
#' @param dge DGEList object with normalization and dispersion already calculated. This should be the output from the "Preprocess_Samba" function.
#' @param method character, indicating which statistical model and test to use for guide-level analysis. Currently, only "QLF" and "LRT" are supported.
#' @param design design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.
#' @param coefficient character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.
#' @param contrast contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).
#' @param GuideMap data.frame that maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.
#' @param file.prefix (optional) file.path and/or character prefix to use for exporting the result tables
#' @return data.frame of the results from the guide-level screen enrichment analysis.
#' @export
Analyze_Samba_Guides <- function(dge, method = 'QLF', design = NULL, coefficient = NULL, contrast = NULL, GuideMap = NULL, file.prefix = NULL){
    # set defaults
    hdr.gene = 'Gene'
    hdr.guide = 'sgRNA'
    fdr.threshold = 0.1

    # check inputs
    if(!xor(is.null(coefficient), is.null(contrast))) warning('Please input a valid coefficient or contrast!')
    if(!(method %in% c('QLF','LRT'))) warning('Please choose either "LRT" or "QLF" as an analysis method!')
    if(is.null(GuideMap)) GuideMap <- dge$GuideMap
    if(is.null(GuideMap)) warning('Please provide a dataframe that maps guides to genes!')

    # get design info
    if(is.null(design)) design <- dge$Design

    # sgRNA-level analysis
    cat('Performing sgRNA-level analysis ...\n')
    if(method == 'QLF'){
        fit <- edgeR::glmQLFit(dge, design = design, robust = T)
        sgRes <- edgeR::glmQLFTest(fit, coef = coefficient, contrast = contrast)
    }
    if(method == 'LRT'){
        fit <- edgeR::glmFit(dge, design = design, robust = T)
        sgRes <- edgeR::glmLRT(fit, coef = coefficient, contrast = contrast)
    }
    cat('Done!\n')

    # Format sg-level results
    sgRes <- edgeR::topTags(sgRes, n = Inf, sort.by = 'PValue')$table
    colnames(sgRes)[3] <- 'Statistic'
    sgRes$sgRNA <- rownames(sgRes)
    sgRes <- merge(GuideMap, sgRes, by = 'sgRNA', all.x = F, all.y = T)
    sgRes$zscore <- scale(sgRes$logFC, center = F, scale = T)[,1]
    sgRes <- sgRes[with(sgRes, order(Statistic,decreasing = T)),]

    # save data
    if(!is.null(file.prefix)){
        saveRDS(fit, file = paste0(file.prefix, '_FittedData.rds'))
        write.table(sgRes, file = paste0(file.prefix, '_GuideLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }

    # return data
    return(sgRes)
}

#' Analyze_Samba_Genes
#'
#' This function performs the gene-level screen enrichment analysis of SAMBA, using guide-level SAMBA results as an input.
#'
#' @param sgRes data.frame of the results from the guide-level screen enrichment analysis. This is the output from the "Analyze_Samba_Guides" function.
#' @param ntc.as.null.dist logical, indicating whether to ese NTCs or all guides to generate a null distribution for gene-level analyses.
#' @param score.method character, indicating whether to use the 'MetaAnalysis' or 'GeneScore' method to perform gene-level analysis.
#' @param file.prefix (optional) file.path and/or character prefix to use for exporting the result tables
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
Analyze_Samba_Genes <- function(sgRes, ntc.as.null.dist = T, score.method = 'MetaAnalysis', file.prefix = NULL){
    # Set defaults
    control.gene <- 'NTC'
    hdr.gene <- 'Gene'
    fdr.threshold <- 0.1
    
    # Check inputs
    if(!(score.method %in% c('GeneScore','MetaAnalysis'))) warning('Please choose either "GeneScore" or "MetaAnalysis" as a score.method!')
    
    ## Filter guide-level data (>= 2sgRNA/gene and not NA)
    sgPerGene <- plyr::ddply(sgRes, 'Gene', nrow)
    keep.genes <- sgPerGene[which(sgPerGene[,2] > 1), 'Gene']
    sgRes <- sgRes[which(sgRes$Gene %in% keep.genes),]
    sgRes <- sgRes[!is.na(sgRes$Gene ),]
    
    ## set data for null dist
    cat('Generating null distribution.\n')
    ifelse(ntc.as.null.dist, nullData <- sgRes[which(sgRes$Gene  == control.gene),], nullData <- sgRes)
    sgRes <- sgRes[which(sgRes$Gene != control.gene),]
    nGeneGuides = sgPerGene[which(sgPerGene$Gene != control.gene),'V1']
    rIndices <- RandomIndexGenerator(sgPerGene = nGeneGuides, nTotalGuides = nrow(nullData))
    nullData <- nullData[unlist(rIndices),]
    i <- lapply(rIndices, length)
    nullData[,hdr.gene] <- paste0('NTC', lapply(1:length(i), function(x){ rep(x,i[[x]]) }) %>% unlist())
    
    # perform enrichment and depletion analyses
    if(score.method == 'GeneScore'){
        cat('Calculating enrichment.\n')
        enriched <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1)
        cat('Calculating depletion.\n')
        depleted <- GeneScoreSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = -1)
        cat('Done!\n')
        output <- merge(enriched,depleted,by = 'Gene')
    }
    if(score.method == 'MetaAnalysis'){
        cat('Calculating enrichment.\n')
        output <- MetaAnalysisSamba(data = sgRes, nullData = nullData, fdr.threshold = fdr.threshold, direction = 1)
        cat('Done!\n')
    }
    output <- output[order(output$pval_pos, decreasing = F),]
    
    # save data
    if(!is.null(file.prefix)){
        write.table(output, file = paste0(file.prefix, '_GeneLevelResults.txt'), sep = '\t', quote = F, row.names = F)
    }
    
    return(output)
}


#' MetaAnalysisSamba
#'
#' This function performs SAMBA gene-level aggregation, using the meta-analysis method. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of guide-level screen enrichment data for the null distribution.
#' @param fdr.threshold numeric between 0 and 1, indicating the FDR to use as a threshold.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
MetaAnalysisSamba <- function(data, nullData, fdr.threshold, direction = 1){
    data$zscore <- direction * data$zscore
    sg.threshold <- quantile(nullData$zscore, 1-fdr.threshold)

    data <- data[with(data, order(zscore, decreasing = T)),]
    bottom <- data[duplicated(data$Gene),]
    bottom <- bottom[duplicated(bottom$Gene),]
    top2 <- setdiff(data$sgRNA, bottom$sgRNA)
    topfdr <- data[which(data$zscore > sg.threshold),'sgRNA']
    sg.keep <- union(topfdr, top2)
    data <- data[which(data$sgRNA %in% sg.keep),]
    fdr.sg <- plyr::ddply(data[which(data$sgRNA %in% topfdr),],'Gene',nrow)
    colnames(fdr.sg) <- c('Gene','FDRGuides')

    ## Gene-level aggregation
    geneRes <- plyr::ddply(data, 'Gene', summarise, score = mean(zscore), fdrGuides = length(Gene), pval = suppressWarnings(metap::sumlog(PValue)$p))
    geneRes$FDR <- p.adjust(geneRes$pval, method = 'fdr')
    geneRes <- merge(geneRes, fdr.sg, by = 'Gene', all.x = T, all.y = F)
    geneRes[which(is.na(geneRes$FDRGuides)),'FDRGuides'] <- 0
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(geneRes) <- c('Gene',paste0('score_',direction),paste0('FilteredGuides_',direction), paste0('pval_',direction),paste0('qval_',direction), paste0('fdrGuides_',direction))

    return(geneRes)
}

#' GeneScoreSamba
#'
#' This function performs SAMBA gene-level aggregation, using the gene score method. Note that this is an internal function.
#'
#' @param data data.frame of the results from the guide-level screen enrichment analysis.
#' @param nullData data.frame of guide-level screen enrichment data for the null distribution.
#' @param fdr.threshold numeric between 0 and 1, indicating the FDR to use as a threshold.
#' @param direction integer of either 1 or -1, indicating whether to perform enrichment or depletion analysis, respectively.
#' @return data.frame of the results from the gene-level screen enrichment analysis.
#' @export
GeneScoreSamba <- function(data, nullData, fdr.threshold, direction = 1){
    data$logFC <- direction * data$logFC
    sg.threshold <- quantile(nullData$logFC, 1-fdr.threshold)
    cat(paste0(c('Enrichment','Depletion')[as.integer(direction == -1)+1],' threshold: ', sg.threshold,'\n'))
    
    # Get Null scores
    nullScores <- WtSumScore(nullData, sg.threshold)

    # Get gene scores
    targetScores <- WtSumScore(data, sg.threshold)

    # Calculate z-score from null distribution
    targetScores$score <- scale(targetScores$score, center = median(nullScores$score), scale = sd(nullScores$score))[,1]
    nullScores$score <- scale(nullScores$score, center = median(nullScores$score), scale = T)[,1]
    
    # Calculate p & q values
    pval <- pnorm((targetScores$score), lower.tail = F)
    qval <- p.adjust(pval, 'fdr')
    
    # Generate output
    output <- data.frame(targetScores, pval = pval, qval = qval)
    ifelse(direction == 1, direction <- 'pos', direction <- 'neg')
    colnames(output)[-1] <- paste0(colnames(output)[-1],'_',direction)
    return(output)
}

#' FilterCountData
#'
#' This function filters DGEList of raw count data, based on a minimum number of total counts.
#'
#' @param dge DGEList object with raw count data.
#' @param group factor, indicating the grouping of all samples. Default is NULL.
#' @param min.guides numeric, indicating the minimum number of counts for each guide across all samples. Default is 0.
#' @return character vector of the guides that passed filtering.
#' @export
FilterCountData <- function(dge, group = NULL, min.guides){
    cat('Filtering guides with low representation across screen samples.\n')
    if(!is.null(group)){
        samples.screen <- which(group != 0)
    } else {
        samples.screen <- colnames(dge$counts)
    }
    keep.exprs = dge$counts[,samples.screen] %>% rowSums(.)
    keep.exprs <- keep.exprs[keep.exprs > min.guides] %>% names()
    cat(paste0('   Removed ',nrow(dge) - length(keep.exprs),' of ',nrow(dge) ,' guides.\n'))
    return(keep.exprs)
}

#' WeighGuides
#'
#' This function calculates a weight for each guide, based on the number of screen samples that have detected counts.
#'
#' @param dge DGEList object with raw count data.
#' @return DGEList object that contains guide weights as a numeric vector.
#' @export
WeighGuides <- function(dge){
    df <- edgeR::cpm(dge, log = F)
    upper_rate_limit <- 1
    min_guide_for_det <- max(quantile(rowMeans(df),.05),1)
    wt <- apply(df,1, function(x) {min(length(x[x > min_guide_for_det])/ncol(df),upper_rate_limit)})
    wt <- ((wt) + .1)
    wt <- 1/(1+exp(-wt))
    wt <- matrix(rep(wt, ncol(dge)), nrow = nrow(dge), byrow = F)
    dimnames(wt) <- dimnames(dge)
    dge$weights = wt
    return(dge)
}


#' RandomIndexGenerator
#'
#' This function creates a list of reproducible, random integers from which to create a null distribution.
#'
#' @param sgPerGene integer vector, representing the number of guides for each gene.
#' @param nTotalGuides integer, indicating the total number of guides in the filtered library.
#' @return list of integer vectors, each of random numbers from which to construct a null data distribution.
#' @export
RandomIndexGenerator <- function(sgPerGene, nTotalGuides){
    # make sure there is between 10,000 and 100,000 iterations
    nIterations <- min(max(length(sgPerGene), 10000), 100000) 
    # make sure there is a min of 2 grna/gene
    sgPerGene[which(sgPerGene < 2)] <- 2
    # if n-genes < n-iterations, bootstrap up to n-iterations
    set.seed(42)
    if(length(sgPerGene) < nIterations) sgPerGene <- dqrng::dqsample(sgPerGene, size = nIterations, replace = T)
    
    set.seed(42)
    randomSeeds <- dqrng::dqsample(1:nIterations, replace = F)
    randData <- lapply(1:nIterations, function(i) {
        set.seed(randomSeeds[i])
        dqrng::dqsample(1:nTotalGuides, size = min(nTotalGuides,sgPerGene[i]), replace = F)
    })
    return(randData)
}


#' WtSumScore
#'
#' This function calculates a single gene score for a given set of guide log-FC values.
#'
#' @param scores numeric vector of the log-FC values for all guides of a single gene.
#' @param sg.threshold numeric, giving the logFC value of the null data that represents the FDR threshold
#' @return numeric score of a gene.
#' @export
WtSumScore <- function(data, fdr.thresh){
    wt.increment <- 1/4
    # sort grna-level results by logfc
    data <- data[with(data, order(Gene, -logFC, decreasing = F)),]
    # get the sum of the top half of grna log-FCs, using quantiles to account for difference in the #s of gRNAs/Gene
    score_tophalf <- data %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::reframe(n.grna = 5, 
                  score_tophalf = sum(quantile(logFC,seq(.5,1,0.125)))/n.grna) 
    # calculate weights for filtered grna
    #data <- data[with(data, order(Gene, -logFC, decreasing = F)),]
    if(nrow(data[which(data$logFC > fdr.thresh),]) > 0){
        score_fdr <- data[which(data$logFC > fdr.thresh),] %>% 
            dplyr::group_by(Gene) %>% 
            dplyr::mutate(wt = 1 + (wt.increment * (1:length(logFC) - 1))) 
	score_fdr[which(score_fdr$wt > 2),'wt'] <- 1
    } else {
        score_fdr <- data %>% 
            dplyr::group_by(Gene) %>% 
            dplyr::mutate(wt = 0) 
    }
    # summarize weighted logFCs
    #score_fdr <- merge(score_fdr, score_tophalf[,c('Gene', 'n.grna')], by = 'Gene')
    score_fdr <- score_fdr %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(n_fdr_guides = length(logFC), score_fdr = sum(logFC * wt)/length(logFC))
    # count the number of grna that pass threshold
    #count_fdr <- plyr::count(data[which(data$logFC > fdr.thresh),'Gene'])
    #colnames(count_fdr) <- c('Gene','n_fdr_guides')
    # merge the score_tophalf, score_fdr, and count_fdr
    scores <- merge(score_fdr, score_tophalf, by = 'Gene', all = T)
    #scores <- merge(scores, count_fdr, by = 'Gene', all = T)
    scores[is.na(scores)] <- 0

    # set max multipliers
    scores$mult_fdr <- scores$n_fdr_guides
    scores$mult_fdr[which(scores$mult_fdr > 5)] <- 5
    scores$mult_fdr <- 0.5 * scores$mult_fdr/max(scores$mult_fdr)
    scores$mult_tophalf <- 1

    scores <- data.frame(Gene = scores$Gene, score_fdr = scores$score_fdr, score_tophalf = scores$score_tophalf,
                         score = (scores$score_fdr * scores$mult_fdr) + 
				 (scores$score_tophalf * scores$mult_tophalf), 
                         n_fdr_guides = scores$n_fdr_guides)

    #scores$score <- scale(scores$score, scale = F, center = T) %>% as.numeric()
    return(scores)
}



#' UpperQuantNorm
#'
#' This function normalizes a count matrix by the 75th percentile of each sample.
#'
#' @param X Count matrix with genes as rows and samples as columns.
#' @return Normalized matrix.
#' @export
UpperQuantNorm <- function(X){
  X<-X+0.1
  upperQ<-apply(X,2,function(y) quantile(y, 0.75))
  f<-upperQ/mean(upperQ) # calculate normalization factor
  res<-scale(X,center=FALSE,scale=f) 
  return(res)
}

#' UpperQuantNorm2Step
#'
#' This function performs a 2-step normalization, first with quantile normalization
#' of samples, then by median normalization of gRNA values.
#'
#' @param X Count matrix with genes as rows and samples as columns.
#' @param f Multiplication factor for nornalized data.
#' @return Normalized matrix.
#' @export
UpperQuantNorm2Step <-function(X, f = 50){
  # per gene normalization by Median: pgQ2
  # X: a matrix of data with the multiplication of factor (f) as:
  # f=100 as a default

  uq.res <- UpperQuantNorm(X) #  perform UQ normalization per sample
  
  ## perform per gene nomralization
  m<-apply(uq.res,1, median)
  
  idx<-which(m==0) # avoid 0 median 
  m[idx]=0.1
  
  si<-m/f # calculate normalization factor
  X1<-scale(t(uq.res),center=FALSE,scale=si)
  res<-t(X1)
  return(res)  
}



                
                
