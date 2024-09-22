#' Calculates CpG-level statistics and identifies stochastic epigenetic mutations (SEMs) from sequencing-based data
#' @param bs A bsseq object
#' @param type Type of methylation values to be extracted from the bsseq object. Either raw or smooth values. For smoothed values, bsseq object needs to have been smoothed prior to running this option.
#' @param minCov Numeric value. The minimum number of reads to overlap a CpG sites by all samples for it to be used for downstream analysis. Note, lowering the number will likely increase the number of CpGs that are tested for SEM identification but at two big costs: 1) increase the computation time as it will be required to analyze more CpGs and 2) likely introduce false positives as less coverage generally equals less certainty of the estimation of methylation at a given CpG site. Default is 10. 
#' @param minSamples Numeric value. The minimum number of samples required to perform the SEM analysis. A larger sample size (e.g., >50) is generally required for an analysis such as this. Default is 50.
#' @param runStats Logical. Whether to run CpG-level statistics (and save them) or not. Default is FALSE.
#' @param nCG Numeric. The number of CpGs used to split data into chunks for easier processing. Default is 100000 CpGs per chunk.
#' @param saveIndSEMs Logical. Whether sample-specific SEMs should be saved. Note, these files may be large in size, depending on how many SEMs are identified per sample. Default is FALSE.
#' @param saveDir Directory where outputs should be saved. Default is current working directory.
#' @param verbose Logical. Whether output of functions (e.g., number of CpGs assess per sample) should be verbose or not. Default is TRUE.
#' @examples
#' findSEMs(bs, minCov = 10, minSamples = 50, runStats = FALSE, cluster = FALSE, nCG = 100000, saveIndSEMs = FALSE, saveDir = getwd(), verbose = TRUE)
#' @import bsseq
#' @import moments
#' @export

findSEMs <- function(bs, type = c("raw", "smooth"), minCov = 10, minSamples = 50, runStats = FALSE, nCG = 100000, saveIndSEMs = FALSE, saveDir = getwd(), verbose = TRUE) {

    # check if sample size is large enough
    if (ncol(bs) < minSamples) {
        stop(paste0("Number of samples is less than "),minSamples,"! Cannot calculate SEMs with such few samples.")
    }

    # filter CpGs based on coverage
    message("[findSEMs]: filtering low-coverage CpGs from further analysis")
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")<minCov) == 0)
    bs <- bs[loci.idx,]

    # get matrix of methylation values
    message("[findSEMs]: getting methylation estimates from CpGs")
    meth.mat <- bsseq::getMeth(bs, type = type)
    rownames(meth.mat) <- paste0(seqnames(bs), ":", start(bs))

    if (runStats == TRUE) {

        cgStats <- function(data) {
            c(q1 = stats::quantile(data, probs = .25, na.rm=TRUE),
            q3 = stats::quantile(data, probs = .75, na.rm=TRUE),
            skewness = moments::skewness(data, na.rm=TRUE),
            sd = stats::sd(data, na.rm=TRUE),
            min = min(data, na.rm=TRUE),
            mean = mean(data, na.rm=TRUE),
            max = max(data, na.rm=TRUE),
            kurtosis = moments::kurtosis(data, na.rm=TRUE),
            iqr = stats::IQR(data),
            range = max(data) - min(data))
        }

        message("[findSEMs]: getting CpG-level statistics now")
        message("[findSEMs]: partitioning data for easier processing")
        nSegments <- round(nrow(meth.mat)/nCG)
        message(paste0("[findSEMs]: splitting methylation data into chunks with ", nCG, " CpGs each"))
        message(paste0("[findSEMs]: a total of ", nSegments, " will be processed"))
        nRows <- 1:nrow(meth.mat)
        x <- split(nRows, factor(sort(rank(nRows)%%nSegments)))
        cpgStats <- c()
        for (ii in 1:length(x)) {
            meth.mat.sub <- meth.mat[x[[ii]],]
            resSub <- apply(meth.mat.sub, MARGIN = 1, FUN = cgStats)
            resSub <- as.data.frame(t(resSub))
            cpgStats <-rbind(cpgStats, resSub)
            if (verbose == TRUE) {
                message(paste0("\tDone with chunk ", ii," of ", length(x)))
            }
        }
 
       message("[findSEMs]: identifying CpG clusters")
        kk <- kmeans(scale(na.omit(cpgStats)), centers = 3)
        cpgStats$cluster[which(!is.na(cpgStats$kurtosis))] <- kk$cluster
 
        message(paste0("[findSEMs]: will save CpG-level statistics at ", saveDir, "cpgStats.rda"))
        save(cpgStats, file = paste0(saveDir, "/cpgStats.rda"))
        rm(nSegments, nRows, x, resSub, meth.mat.sub)
        gc()

    }

    # workhorse . . . time to find SEMs
    message("[findSEMs]: finding SEMs now . . .")
    SEMFinder <- function(data) {
        hyperSEM <- 0
        hypoSEM <- 0
        iqr <- stats::IQR(data)
        q3 <- stats::quantile(data, probs = .75, na.rm=TRUE)
        q1 <- stats::quantile(data, probs = .25, na.rm=TRUE)
        extremeHyper <-  q3 + (3 * iqr)
        extremeHypo <- q1 - (3 * iqr)
        if (iqr > 0) {
            ix.hyper <- which(data > extremeHyper)
            ix.hypo <- which(data < extremeHypo)
            if (length(ix.hyper) > 0) { 
                totalSEMs[ix.hyper, 1] <-  totalSEMs[ix.hyper, 1] + 1
            }
            if (length(ix.hypo) > 0) { 
                totalSEMs[ix.hypo, 2] <-  totalSEMs[ix.hypo, 2] + 1
            }
        }
        return(totalSEMs)
    }

    nSegments <- round(nrow(meth.mat)/nCG)
    message(paste0("[findSEMs]: splitting methylation data into chunks with ", nCG, " CpGs each"))
    message(paste0("[findSEMs]: a total of ", nSegments, " chunks will be processed"))
    nRows <- 1:nrow(meth.mat)
    x <- split(nRows, factor(sort(rank(nRows)%%nSegments)))
    semTotals <- c()
    for (ii in 1:length(x)) {
         totalSEMs <- data.frame(matrix(0, nrow = ncol(meth.mat), ncol = 2))
        meth.mat.sub <- meth.mat[x[[ii]],]
        resSub <- apply(meth.mat.sub, MARGIN = 1, FUN = SEMFinder)
        semTotals <- c(semTotals, resSub)
        if (verbose == TRUE) {
            message(paste0("\tDone with chunk ", ii," of ", length(x)))
        }
    }
    # clean ‘er up
    message("[findSEMs]: bringing it all back home")
    message("\tSummarizing hyper-SEMs")
    hyperCount <- semTotals[[1]][1]
    for (i in 2:length(semTotals)) {
    hyperCount <- hyperCount + semTotals[[i]][1]
    if(abs(i) %% 1000000 == 0) {
        message(paste0("\t\t Done with ", i, " CpGs . . ."))
        }
   }
    hypoCount <- semTotals[[1]][2]
    message("\tSummarizing hypo-SEMs")
    for (i in 2:length(semTotals)) {
        hypoCount <- hypoCount + semTotals[[i]][2]
    if(abs(i) %% 1000000 == 0) {
        message(paste0("\t\t Done with ", i, " CpGs . . ."))
        }
    }
    sampleSEMs <- cbind(hyperCount, hypoCount)
    colnames(sampleSEMs) <- c("Hyper_SEMs" ,"Hypo_SEMs")
    sampleSEMs$Total_SEMs <- sampleSEMs$Hyper_SEMs + sampleSEMs$Hypo_SEMs
    sampleSEMs$log10SEMs <- log10(sampleSEMs$Total_SEMs)
    return(sampleSEMs)
}



#' Identify low methylated regions (LMRs) (and unmethylated regions (UMRs)) for each sample from a bsseq object
#' @param bs A bsseq object
#' @param minCov Minimum coverage by all samples across a CpG to be considered for further analysis. Default is 10
#' @param methCutoff Methylation value cutoff (from 0 to 1) used to look for LMRs/UMRs. Default is 0.5 (50%)
#' @param numCG Minimum number of CpGs used to find LMRs/UMRs. Default is 4
#' @param chrSelect The chromosome used to train the partially methylated domains (PMDs) hidden Markov models. A shorter autosome is generally recommended. Default is chr22
#' @param umrs Logical. Whether or not to return UMRs. Default is FALSE.
#' @param nCores Number of cores to use for PMD-HMM training and LMR/UMR identification. Default is 1
#' @examples
#' findLMRs(bs, methCutoff = 0.5, numCG = 4, chrSelect = "chr22", umrs = FALSE, nCores = 1)
#' @import bsseq
#' @import GenomicRanges
#' @import MethylSeekR
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export

findLMRs <- function(bs, minCov = 10, methCutoff = 0.5, numCG = 4, chrSelect = "chr22", umrs = FALSE, nCores = 1) {

    message("[findLMRs]: filtering low-coverage CpGs from further analysis")
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type = "Cov") < minCov) == 0)
    bs <- bs[loci.idx,]

    chrLengths <- seqlengths(Hsapiens)

    results <- list()

    for (i in 1:ncol(bs)) {

        message(paste0("Starting LMR identification for sample ", i))

        oneSample <- GenomicRanges::GRanges(seqnames(bs[,i]), IRanges(start(bs[,i])))
        mcols(oneSample)$T <- bsseq::getCoverage(bs[,i], type = "Cov")
        mcols(oneSample)$M <- bsseq::getCoverage(bs[,i], type = "M")

        pmds.gr <- MethylSeekR::segmentPMDs(m = oneSample, chr.sel = chrSelect, seqLengths = chrLengths, num.cores = nCores)

        sampleLMRs.gr <- MethylSeekR::segmentUMRsLMRs(m = oneSample, meth.cutoff = methCutoff, nCpG.cutoff = numCG, PMDs = pmds.gr, num.cores = nCores, myGenomeSeq = Hsapiens, seqLengths = chrLengths)

        if (umrs == FALSE) {
            sampleLMRs.gr <- sampleLMRs.gr[which(sampleLMRs.gr$type == "LMR"),]
        }

        results[[i]] <- sampleLMRs.gr

    }

    return(results)

}



#' Identify imprinted control regions (ICRs) from a bsseq object
#' @param bs A bsseq object
#' @param minCov Minimum coverage of a CpG from a sample to be considered for further analysis. Default is 1
#' @param minSample Minimum portion of samples to meet minCov cutoff to be considered for further analysis. Default is 0.5 (50%)
#' @param numCG Minimum number of adjacent CpGs used to be tested in the sliding window to identify ICRs. Default is 5
#' @param tissueType Type of tissue of samples. One of "somatic" or "gamete". Default is "somatic".
#' @param lowerLimitSomatic Lower limit of methylation percentage to use to determine hemimethylation, if tissue is somatic. Default is 0.35 (35%). Note, this cutoff is recommended for somatic tissues.
#' @param upperLimitSomatic Lower limit of methylation percentage to use to determine hemimethylation, if tissue is somatic. Default is 0.65 (65%). Note, this cutoff is recommended for somatic tissues. 
#' @param lowerLimitGamete Lower limit of methylation percentage to use to determine hemimethylation, if tissue is gametic. Default is 0.35 (35%). Note, this cutoff is recommended for somatic tissues.
#' @param upperLimitGamete Lower limit of methylation percentage to use to determine hemimethylation, if tissue is gametic. Default is 0.65 (65%). Note, this cutoff is recommended for somatic tissues. 
#' @examples
#' findICRs(bs, minCov = 1, minSample = 0.5, numCG = 5, lowerLimit = 0.35, upperLimit = 0.65)
#' @import bsseq
#' @import GenomicRanges
#' @export

findICRs <- function(bs, minCov = 1, minSample = 0.5, numCG = 5, tissueType = NULL, lowerLimitSomatic = 0.35, upperLimitSomatic = 0.65, lowerLimitGamete = 0.1, upperLimitGamete = 0.8) {

    if (is.null(tissueType)) {
        tissueType <- "somatic"
    }
    if ((tissueType != "somatic") & (tissueType != "gamete")) {
        stop("tissueType must be one of gamete or somatic")
    }
    if (tissueType == "somatic") {
        lowerLimit <- lowerLimitSomatic
        upperLimit <- upperLimitSomatic
    }
    if (tissueType == "gamete") {
        lowerLimit <- lowerLimitGamete
        upperLimit <- upperLimitGamete
    }
    message(paste0("[findICRs]: will look for ICRs from ", tissueType, " using limits of [", lowerLimit,", ", upperLimit, "]"))
    message("[findICRs]: filtering low-coverage CpGs from further analysis")

    ix <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") < minCov) < (ncol(bs) * minSample))
    bs <- bs[ix,]

    message("[findICRs]: getting methylation percentages now")
    meth.mat <- bsseq::getMeth(bs, type = "raw")
    rownames(meth.mat) <- paste0(seqnames(bs), ":", start(bs))

    message("[findICRs]: identifying ICRs now")
    candidate_icrs <- c()
    currentChr <- ""
    for (i in 1:(nrow(meth.mat)-(numCG - 1))) {

        ii <- i + (numCG - 1)

        meth.sub <- meth.mat[i:ii,]

        meth.sub2 <- as.data.frame(meth.sub)
        icr.chr <- as.data.frame(t(as.data.frame(strsplit(rownames(meth.sub2)[1], ":"))))
        if (icr.chr[1,1] != currentChr) {
            currentChr <- icr.chr[1,1]
            message(paste0("\tWorking on ", currentChr))
        }

        if (tissueType == "somatic") {
            if (min(meth.sub, na.rm=TRUE) >= lowerLimit & max(meth.sub, na.rm=TRUE) <= upperLimit) {
                meth.sub <- as.data.frame(meth.sub)
                icr.start <- rownames(meth.sub[1,])
                icr.end <- rownames(meth.sub[numCG,])
                icr.start <- as.data.frame(t(as.data.frame(strsplit(icr.start, ":"))))
                icr.end <- as.data.frame(t(as.data.frame(strsplit(icr.end, ":"))))
                icr.x <- cbind(icr.start, icr.end)
                colnames(icr.x) <- c("chr", "start", "chr2", "end")
                rownames(icr.x) <- "ICR"
                if (icr.x$chr == icr.x$chr2) {
                    candidate_icrs <- rbind(candidate_icrs, icr.x)
                }
            }
        }
        if (tissueType == "gamete") {
            if (max(meth.sub, na.rm=TRUE) <= lowerLimit | min(meth.sub, na.rm=TRUE) >= upperLimit) {
                meth.sub <- as.data.frame(meth.sub)
                icr.start <- rownames(meth.sub[1,])
                icr.end <- rownames(meth.sub[numCG,])
                icr.start <- as.data.frame(t(as.data.frame(strsplit(icr.start, ":"))))
                icr.end <- as.data.frame(t(as.data.frame(strsplit(icr.end, ":"))))
                icr.x <- cbind(icr.start, icr.end)
                colnames(icr.x) <- c("chr", "start", "chr2", "end")
                rownames(icr.x) <- "ICR"
                if (icr.x$chr == icr.x$chr2) {
                    candidate_icrs <- rbind(candidate_icrs, icr.x)
                }
            }
        }
    }

    message("[findICRs]: cleaning up now")
    if (is.null(candidate_icrs) == TRUE) {
        message("[findICRs]: no ICRs were identified. Bummer . . .")
    } else {    
        candidate_icrs <- as.data.frame(candidate_icrs)
        colnames(candidate_icrs) <- c("chr", "start", "chr2", "end")
        candidate_icrs.gr <- with(candidate_icrs, GRanges(chr, IRanges(as.numeric(start),as.numeric(end))))
        candidate_icrs <- as.data.frame(reduce(candidate_icrs.gr))
        return(candidate_icrs)
    }
}

