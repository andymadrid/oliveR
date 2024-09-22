## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  ### Install required packages
#  install.packages("devtools")
#  library(devtools)
#  
#  ### install the oliveR package
#  install_github("andymadrid/oliveR")
#  
#  

## -----------------------------------------------------------------------------
# load in the oliveR pacakge
suppressPackageStartupMessages(library(oliveR))



## -----------------------------------------------------------------------------
# get toy example
data(bs)



## -----------------------------------------------------------------------------
# summary of bsseq object
bs

# summary of phenotypic data associated with bsseq object
pData(bs)



## -----------------------------------------------------------------------------
# guess sample sex from bsseq object
guessed_sex <- guessSex(bs)

# look at the output
guessed_sex



## -----------------------------------------------------------------------------
# predict age from bsseq object
predicted_age <- madrid(bs)

# look at the output
predicted_age



## -----------------------------------------------------------------------------
# assess the correlation between actual and estimated age
corMADRID <- cor(predicted_age$MADRID_Age, pData(bs)$Age)
corMADRID


## -----------------------------------------------------------------------------
# test for differences (accelerated/decelerated) in actual and estimated age
predicted_age$difference <- predicted_age$MADRID_Age - pData(bs)$Age

t.test(predicted_age$difference ~ factor(pData(bs)$Diagnosis))



## -----------------------------------------------------------------------------
# test for differences in the resduals between actual and estimated age
predicted_age$residuals <- resid(lm(predicted_age$MADRID_Age ~ pData(bs)$Age))
t.test(predicted_age$residuals ~ factor(pData(bs)$Diagnosis))



## ----eval = FALSE-------------------------------------------------------------
#  # estimate telomere lengths
#  tls <- estimateTL(bs)
#  
#  # look at the output
#  tls
#  

## ----eval=FALSE---------------------------------------------------------------
#  # estimate cell-type proportions from blood samples
#  cellProps <- estimateCellProps(bs)
#  
#  # check the estimates
#  head(cellProps)
#  

## -----------------------------------------------------------------------------
# lift coordinates from hg38 to hg19
bsLifted <- liftBuild(bs, current = "hg38", new = "hg19")

# check the lifted coordinates
bsLifted


## -----------------------------------------------------------------------------
# load in DMR GRanges object
data(dmrs)

# take a look at the DMRs
dmrs

# overlap with the imprintome
dmrsImprinted <- imprintMe(dmrs)

# look at the output
dmrsImprinted


## ----eval=FALSE---------------------------------------------------------------
#  # identify ICRs in toy dataset
#  icrs <- findICRs(bs)
#  
#  # take a look at the ICRs
#  icrs
#  

## -----------------------------------------------------------------------------
# detect SEMs
#sampleSEMs <- findSEMs(bs, minSamples = 10) # setting samples to 10 for demo, but >50 is recommended for actual use


# check the output
sampleSEMs

# check to see if there are any differences between groups
t.test(sampleSEMs$log10SEMs ~ factor(pData(bs)$Diagnosis))


## ----eval=FALSE---------------------------------------------------------------
#  # Identify LMRs
#  lmrs <- findLMRs(bs)
#  
#  # take a look at LMRs
#  lmrs
#  

## -----------------------------------------------------------------------------
# guess smoking status
smoking <- guessSmoking(bs)

# look at the results
smoking


