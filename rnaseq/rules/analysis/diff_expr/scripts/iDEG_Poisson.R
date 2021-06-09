## a function to detect DEG when assuming RNA-Seq data follows Poisson distribution
## Qike Li 5/23/2016


####### define a function for finding DEGs under negative binomial distribution and non-constant overdispersion assumptions 
iDEG.Pois <- function (dataTest, data1, ind0,
                        normalization, df, nulltype, pct, plot){ #dataTest needs to be count matrix, baseline needs to be baseline, case is the case

## modified TPM normalizatioin
#if(normalization) dataTest <- apply(dataTest,2,function(x) x/sum(x)*10^6) # NOTE: the number 10^6 does affect the result, but not too much, may consider to multiply the mean of the two library size, instead of a fixed number 10^6
if(normalization) dataTest[,'case'] <- dataTest[,'case'] * (sum(dataTest[,'baseline'])/sum(dataTest[,'case']))

    ## Varaice stabilizing transformation for Poisson data_
    data_vst <- f.vst.pois(dataTest)
    ## use median absolute deviation to estimate the variance
    sd_hat <- stats::mad(data_vst[,2]-data_vst[,1])
    ## conduct Z test to give every gene a p-value
    zz <- (data_vst[,2]-data_vst[,1])/sd_hat
    ## hist(zz,breaks=110)
    local_fdr <- locfdr::locfdr(zz, df = df, nulltype = nulltype, plot = plot, pct = pct)$fdr
    data1 <- as.data.frame(data1)
    data1[!ind0,'local_fdr'] <- local_fdr
    data1[ind0,'local_fdr']  <- 1 # we assign local_fdr 1 for the genes whose expression levels are zero in both conditions
    data1[!ind0,'statistics'] <- zz
    data1[ind0,'statistics'] <- 0
    ls_tmp <- list(result = data1) 
    return(ls_tmp)
}

f.vst.pois <- function(x) sqrt(x) + sqrt(x+1)
