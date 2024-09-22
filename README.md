# oliveR: an R package of DNA methylation microarray functionalities for sequencing-based data

## Welcome
Hey, thanks for stopping on by! This is my first R package, so bear with me here.

The package is pretty straightforward; all you really need to start is a bsseq object.

As DNA methylation microarrays have been around for awhile now, several functionalities and predictive models have been generated for microarray datasets. However, as sequencing-based DNA methylation datasets (e.g., whole-genome bisulfite sequencing) have not been around as long - mainly due to cost restraints - these same functionalities that were developed for microarray data have not been implemented in sequencing datasets. As such, this R packages seeks to bridge that gap by allowing those array-based utilities to be functional with sequencing data. Currently, starting with just a bsseq object, you can estimate sample ages, identify stochastic epigenetic mutations (SEMs), find lowly-methylated regions (LMRs), find putative imprinted control regions (ICRs), guess the sex of samples, estimate telomere length, and more!

## Vignette

Check out the tutorial of the general functionality of the package at:

http://htmlpreview.github.io/?https://raw.githubusercontent.com/andymadrid/oliveR/main/vignettes/introduction.html

## Sequencing-based data to array-based clocks

Additional functionality has been added to the madrid() function to now allow it to predict age from several array-based clocks, including Horvath's (original), Hannum's, Levine's (PhenoAge), Horvath's (Skin-Blood), Lin's, and DunedinPACE. Essentially, it takes a matrix of methylation values, formats them a bit, filters to only those CpGs on the arrays (27k, 450k, and EPICv1), adds CpG names from the arrays, then leans heavily on R packages wateRmelon and DunedinPACE to estimate age from all those clocks.

## Imputation

Your dataset may very well have some missing CpGs for either the sequencing-based clock and/or the array-based clocks, which is totally okay. It's to be expected, especially if sequencing depth was on the lower end. To get around this, an imputation step has been implemented in the predictAge() function which calls the mice() function to impute missing data, using a random forest method. In my testing I found that - at least for the sequencing-based clock - mice() run with random forest then averaged across imputed datasets yielded better estimates than, say, impute.knn() did, which is the method commonly used for most array-based clocks. An advanced user can feel free to edit the code as they see fit to change the imputation method. The currently employed method just happened to be the best that I found in my testing, but may not be robust enough for all datasets. However, an astute observer will notice that the imputation step for the array-based clocks actually uses impute.knn(). I found that it worked quite well for the array clocks (who would have guessed?) and worked pretty fast.  

## Cell-type proportion estimates

I've also included a function (estimateCellCounts()) that takes your sequencing data and estimates cell-type proportions. To do this, it formats and filters data from a bsseq object then leans heavily on the R package meffil to do the calculation. There are other packages out there that do similar things (e.g., methylcc), but I found - through my own testing, using in-house datasets - that this method produced results that were more correlative with array data; correlations between
estimates from samples run on the EPIC array and the same samples that were sequenced were ~0.95.

## Stochastic epigenetic mutations (SEMs)

While the main function of this package is to estimate age from sequencing-based data, others have recently been moving to detecting stochastic epigenetic mutations (SEMs) in DNA methylation data. So, I've implemented a function to do just that here. It calculates the interquartile range (IQR) for each CpG, then goes sample by sample, CpG by CpG to identify SEMs. For more information on SEMs, please refer to Gentilini et al (2015) Aging and Markov et al (2024) GeroScience. Please note that the function is slow and can be computationally expensive. I recommend larger sample sizes (>50) and higher coverage (>10) to reduce the number of CpGs analyzed and to only keep those with substantial evidence of their estimated methylation levels.

## Word of warning

### Coverage
The biggest issue that I have run into while developing this package has been coverage. For reliable, robust results, higher coverage from sequencing data is needed. While testing, I found several datasets that had ~1x genome-wide coverage. As such, estimates for age, sex, LMRs, SEMs either would highly inaccurate, or were unable to be computed. This is one of the major drawbacks of sequencing data, relative to the array where most samples have useable data on just about every CpG tested on the array. Several imputation methods are utilized in the functions throuhgout the package, but those can only do so much. So, just be wary of that fact. As with most sequencing-based analyses: higher coverage gives better, more robust, and more accurate/reliable results.

### Genome reference
Samples used to test and validate this clock were aligned to the human genome (hg38) using the UCSC chromsome naming scheme (e.g,. chr1, chr2, chr3, etc). If you happen to have aligned your samples to, say, hg19, those coordinates will have to be lifted over in order to properly get filtered/selected for the clock to work. I added in a function liftBuild() that does just that.  

## What's in a name?
You may very well be sitting there (though let's be honest, you're not) and wondering, "what does the name oliveR have to do with DNA methylation?" The answer is a resounding . . . nothing! My dog's name is Olive, and this is an R package. So, through the transitive property, I present you with the R package oliveR! 

## Questions, concerns, collaborations
Do you have your own sequencing-based methylation data that you want to build your own clock with? If so, feel free to email me at at madrid2[at]wisc.edu and we can work it out, together. I am always happy to collaborate and/or help! Also, if there's some functinoality that you think would be of great value to add, you can also let me know and I can work on implementing that, as well.

See you space cowboy...
- AM
