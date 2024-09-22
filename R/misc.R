#' Overlaps DMRs with known imprinted regions
#' @param dmrs A GenomicRanges object containing the chromosome, start, and end of regions of interest (e.g., differentially methylated regions (DMRs)), with coordinates corresponding to the hg38 version of the human genome
#' @examples
#' imprintMe(dmrs)
#' @import GenomicRanges
#' @export

imprintMe <- function(dmrs) {

    # check that the object is a GenomicRanges object
    if (intersect("GRanges", class(dmrs)[1]) == "GRanges") {        
        message("[imprintMe]: Looking for overlaps with imprinted regions now")
    } else {
        stop("[imprintMe]: Your object does not seem to be a GRanges object. Please double check your input object.")
    }

    # load in known imprinted regions
    data(imprintome)
    
    # look for overlaps
    overlaps <- suppressWarnings(as.data.frame(GenomicRanges::findOverlaps(dmrs, imprintome)))

    # filter for the overlaps, then return object
    dmrs.sub <- dmrs[overlaps$queryHits, ]
    imprintome.sub <- imprintome[overlaps$subjectHits, ]    
    mcols(dmrs.sub) <- cbind(mcols(dmrs.sub), as.data.frame(imprintome.sub))
    return(dmrs.sub)
}

#' Inverse transformation of predicted age
#' @param tAge Transformed age from prediction model
#' @param adult.age Adult age for transformation. Default is 20
#' @examples
#' inverse.transform(bs, adult.age = 20)
#' @export
inverse.transform <- function(tAge, adult.age=20) {
        if (tAge < 0) {
                iAge <- (exp(tAge + log(adult.age + 1)))-1
        }
        else {
                iAge <- (tAge*(adult.age+1))+(adult.age)
        }
        return(iAge)
}

#' Transformation of actual age
#' @param age Age of sample(s)
#' @param adult.age Adult age for transformation. Default is 20
#' @examples
#' transform_age(age, adult.age = 20)
#' @export
transform_age <- function(age, adult.age=20) {
        if (age <= adult.age) {
                tAge <- log(age+1)-log(adult.age+1)
        }
        else {
                tAge <- (age-adult.age)/(adult.age+1)
        }
        return(tAge)
}


#' Convert coordinates from hg19 (or hg38) to hg38 (or hg19)
#' @param bs A bsseq object
#' @param current The current human genome build. Can be one of hg19 or hg38 
#' @param new The human genome build that you would like to convert to. Can be one of hg19 or hg38
#' @examples
#' liftBuild(bs, current = "hg38", new = "hg19")
#' @import AnnotationHub
#' @import rtracklayer
#' @import MatrixGenerics
#' @import GenomicRanges
#' @import dplyr
#' @export

liftBuild <- function(bs, current = c("hg19", "hg38"), new = c("hg19", "hg38")) {

    # check the current and new version of the genome you want
    if (current == "hg38" && new == "hg19") {
        direction <- "AH14108"
    }    else if (current == "hg19" && new == "hg38") {
        direction <- "AH14150"
    }  else {
        stop("Currently only can go from hg19 to hg38 (or the reverse). Consider what builds you’re trying to lift")
    }
    message(paste0("[liftBuild]: Will lift your bsseq object from ", current," to ", new))
    # hub <- AnnotationHub()
    # chains <- query(hub, c(current, new, "chainfile"))
    # AH14108 if hg38 to hg19
    # AH14150 if hg19 to hg38
    mcols(bs)$cpgs <- 1:length(bs)
    hgCurrent <- rowRanges(bs)
    hgCurrent$cpgs <- 1:length(hgCurrent)
    hgNew <- hgCurrent %>% rtracklayer::liftOver(AnnotationHub::AnnotationHub()[[direction]]) %>% unlist()
    bsLifted <- bs[which(hgCurrent$cpgs %in% hgNew$cpgs)]
    rowRanges(bsLifted) <- hgNew
    genome(bsLifted) <- new
    lost <- length(bs) - length(bsLifted)
    message(paste0("[liftBuild]: Lifting over resulted in a loss of ", lost," CpGs"))
    return(bsLifted)
}

