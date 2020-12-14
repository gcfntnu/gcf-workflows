## source('~/Dropbox/Qike/adaptive_cutoff/functions/wml source code.R')

####### define a function for finding DEGs under negative binomial distribution and non-constant overdispersion assumptions 
iDEG.NB <- function (dataTest, data1, ind0 ,
                     ## cutOff,
                     normalization, numBin,rankBaseline, 
                     estSize,
                     estBaseline, spar,
                     constDisp, df, nulltype, pct, plot){ #data.test needs to be count matrix, baseline needs to be baseline, case is the case
    if(normalization){
    cds <- edgeR::DGEList(dataTest)
    cds <- edgeR::estimateCommonDisp(cds)
    dataTest <- cds$pseudo.counts    
        }
### order all genes based on expression mean
    if(rankBaseline){      #if rank based on baseline transcriptome
        ## 
        ## bin genes according to percentiles based on the baseline transcriptome
        ind_bins <- split(1:dim(dataTest)[1],cut(dataTest[,'baseline'],breaks = unique(round(stats::quantile(dataTest[,'baseline'],seq(1/numBin,1, by = 1/numBin)),5)),includeLowest = T))
    }else if (!rankBaseline){             #if ran based on the average of the transcritpomes
###########################################################################
        ## note that the differentially expressed genes may have their mean expression shifted up or down
        expMeanHat <- rowMeans(dataTest)    #mean of the two expression levels for each gene
        ## bin genes according to percentiles based on the average of the two transcriptomes
        ind_bins <- split(1:dim(dataTest)[1],cut(expMeanHat,breaks = unique(stats::quantile(expMeanHat,seq(1/numBin,1, by = 1/numBin))),includeLowest = T))
###########################################################################
    }
    
    ## apply to every window the function that estimates the parameters for each window 
    par_bin <- sapply(ind_bins,f.window.mad, dataTest,estBaseline) 

    ## calculate variance - mean for each window
    xTmp <- pmax(par_bin['bin_var',]- par_bin['bin_mean',],0)
    ## calcuate mean^2 for each window
    yTmp <- par_bin['bin_mean',]^2
    
    ## ==== estiamte Size/dispersion parameter 
    if(constDisp){ # if we assume overdispersion is the same across all windows
        ## estimate the constant overdispersion by fitting a linear regression model withou intercept, since variance-mean = mean^2*dispersion
        if(estSize){                    # if estimate size
            sizeHat <- stats::lm(yTmp~0+xTmp)$coefficients
        }else{                          # if estimate dispersion
            dispHat <- stats::lm(xTmp~0+yTmp)$coefficients
        }
    }else{ # if we assume overdispersion vary across windows but similar when expression levels are similar
        if(estSize){                    # if estimate size
            fit_s <- stats::smooth.spline(par_bin['bin_mean',],yTmp/xTmp, spar = spar)
        } else{                     #if estiate dispersion
            fit_s <- stats::smooth.spline(par_bin['bin_mean',],xTmp/yTmp, spar = spar)
        }
        ## size.hat <- predict(fit.s, data.test[,'baseline'])$y
        fitted_y <- stats::predict(fit_s, par_bin['bin_mean',])$y
        fitted_y <- pmax(fitted_y, min(fitted_y[fitted_y>0]))
        if(length(ind_bins) != length(fitted_y)) stop('check alpha function fitting')
        tmpHat <- numeric(length = dim(dataTest)[1])
        for( i in 1:length(ind_bins) ){
            tmpHat[ ind_bins[[i]] ] <- fitted_y[i]
        }
        tmpHat[tmpHat == 0] <- fitted_y[1]
        if(length(tmpHat)!=dim(dataTest)[1]) message('check the number of size estimates')
        if(estSize){
            sizeHat <- tmpHat
        }else{
            dispHat <- tmpHat
        }         
    }
    
    ## Varaice stabilizing transformation for negative binomial data. size is estimated above--size.hat
    if(!estSize) sizeHat <- 1/dispHat
    if(min(sizeHat)<=1.5) {
        data_vst <- cbind(f.nb.vst.2(dataTest[,'baseline'], r = sizeHat),
                          f.nb.vst.2(dataTest[,'case'], r = sizeHat))
    } else {
    data_vst <- cbind(f.nb.vst.1(dataTest[,'baseline'], r = sizeHat),
                     f.nb.vst.1(dataTest[,'case'], r = sizeHat))
    }
    ## use median absolute deviation to estimate the variance
    ## sd_hat <-  mad(data_vst)

    ## compute local fdr
    zz <- (data_vst[,2]-data_vst[,1])/stats::mad(data_vst[,2]-data_vst[,1])
    local_fdr <- locfdr::locfdr(zz, plot = plot, df = df, pct = pct, nulltype = nulltype)$fdr
    data1 <- as.data.frame(data1)
    ## data1[!ind0,'p_value'] <- p_diff
    data1[!ind0,'local_fdr'] <- local_fdr
    data1[ind0,'local_fdr']  <- 1 # we assign local_fdr 1 for the genes whose expression levels are zero in both conditions
    data1[!ind0,'statistics'] <- zz
    data1[ind0,'statistics'] <- 0
    ## DE_status <- rep(FALSE, dim(data1)[1])
    ## DE_status[!ind0] <- data1[!ind0,'local_fdr']<=cutOff
    if(constDisp){
        ls <- list(result = data1, sizeHat = sizeHat, zeroGenes = ind0)
    }else if(!constDisp){
        ls <- list(result = data1, sizeHat = sizeHat, zeroGenes = ind0, alphaFunction = fit_s)
    }
        
    return(ls)
}

####
f.window.mad <- function(index, data_par_est,estBaseline){ # a function to calculate the variance and mean given a group of negative binomial data; x is the indexes of a vector of Negative Binomial gene paired data belonging to the same window
    data_par_est <- data_par_est[index,]
    if(estBaseline){  # if we estimate the parameters for each window based on basline 
        data_tmp <- data_par_est[,'baseline']
        var_bin <- stats::mad(as.numeric(data_tmp))^2
    }else if (!estBaseline){  # if we estimate the parameters for each window based on both transcriptomes
        data_diff <- data_par_est[,'case'] - data_par_est[,'baseline'] # taking the difference of the paired expression value for each gene
        ## ######################################################################################################################
        ## note that here using MAD to estimate standard deviation is not correct. when we use MAD to estimate the standard deviation of Negative Binomial distribution, the constant that need to be multipled to mad to get sd is different from the constant used for Normal distribution
        var_bin <- stats::mad(as.numeric(data_diff))^2/2 # estimate the variance of x, assuming they have the same gene expression mean; data1 needs to be all the differences of every gene's expression values between two conditions. mad(as.numeric(data1[index]))^2 is divided by 2 because: var(X_1 - X_2) = var(X_1) + var(X_2) when X_1 and X_2 are independent 
    }
    ## ######################################################################################################################
    mean_bin <- stats::median(as.numeric(unlist(data_par_est))) # taking the avereage of all values in x as the estimated expression mean of this window }
    pars_oneWin <- c(bin_var = var_bin, bin_mean = mean_bin) # store the variance and mean of x in a data frame
    return(pars_oneWin)
}

f.nb.vst.2 <- function(x,r)  sqrt(r)*asinh(sqrt(x/r))# VST for Negative Binomial distribution  r is overdispersion parameter of a Negative Binomial distribution

f.nb.vst.1 <- function(x,r)  sqrt(r)*asinh(sqrt(x/r)) +
                               sqrt(r - 1)*asinh(sqrt((x+.75)/(r-1.5))) # VST for Negative Binomial distribution  r is overdispersion parameter of a Negative Binomial distribution

## f.window.par.wml <- function(index, data.par.est ){ # a function to calculate the variance and mean given a group of negative binomial data; index is a vector of indexes for a vector of Negative Binomial gene paired data belonging to the same window; data.par.est is the data used for calculating the mean of the genes in one bin
##                                         #data.diff = data.par.est[,2] - data.par.est[,1] # taking the difference of the paired expression value for each gene
## ########################################################################################################################
##     ## note that here using MAD to estimate standard deviation is not correct. when we use MAD to estimate the standard deviation of Negative Binomial distribution, the constant that need to be multiplied to mad to get sd is different from the constant used for Normal distribution
##                                         #var.bin = mad(diff.trim[index])^2 # estimate the variance of x, assuming they have the same gene expression mean 
##                                         #var.bin = mad(as.numeric(data.diff[index]))^2/2 # estimate the variance of x, assuming they have the same gene expression mean; data1 needs to be all the differences of every gene's expression values between two conditions. mad(as.numeric(data1[index]))^2 is divided by 2 because: var(X_1 - X_2) = var(X_1) + var(X_2) when X_1 and X_2 are independent 
## ########################################################################################################################
##                                         #mean.bin = mean(rowMeans(data.trim[index,])) # taking the avereage of all values in x as the estimated expression mean of this window
##     wml.bin.res.tmp <- wml(as.numeric(data.par.est[index,]))
##     mean.bin <- wml.bin.res.tmp$Estimates['m','WML']
##     alpha.bin <- wml.bin.res.tmp$Estimates['alpha','WML']
##     pars.oneWin <- c(alpha = alpha.bin, bin_mean = mean.bin) # store the variance and mean of x in a data frame
##     return(pars.oneWin)
## }

f.window.par.MLE <- function(index, dat.baseline ){ # a function to calculate the variance and mean given a group of negative binomial data; index is a vector of indexes for a vector of Negative Binomial gene paired data belonging to the same window; data.par.est is the data used for calculating the mean of the genes in one bin
    
    ## index <- ind.bins[[1]]
    ## extract the data for one bin
    dat.bin <- dat.baseline[index]
    ## estimate MLEs from the data in this bin
    ## mle.tmp <- fitdistr(dat.bin, densfun = "negative binomial")$estimate
    m <- mean(dat.bin)
    mle.tmp <- c(size = m^2/(stats::var(dat.bin) -m), mu = m)
    pars.oneWin <- c(dispersion = mle.tmp['size'], bin_mean = mle.tmp['mu']) # store the variance and mean of x in a data frame
    return(pars.oneWin)
}
