#' Identification of \strong{i}ndividualized \strong{D}ifferentially
#' \strong{E}xpressed \strong{G}enes (iDEG).
#'
#' Identify differentionally expressed genes between two conditions,
#' and only one transcriptome is collected for each condition.  
#' 
#' @param baseline a vector of gene expression levels of the baseline
#' transcriptome (e.g., healthy tissue)
#' @param case a vector of gene expression levels of the case transcriptome (e.g., tumor tissue)
#' @param normalization a logical variable indicating if normalization has been done
#' @param dataDistribution the distribuitonal assumption of the RNA-Seq data under analysis. Possible values are 'Poisson' and 'NB'. Default is NB--negative binomial.
#' @param numBin number of bins used to group all genes into. Default is 100.
##' @param rankBaseline if True, iDEG groups all genes based on the gene expression levels of the baseline transcriptome. If False, iDEG group all genes based on the gene expression levels of the average of baseline and case transcriptomes.
#' @param constDisp if True, iDEG assumes the dispersion is a count across all genes. If False, iDEG assume dispersion is a smooth fucntion os expression mean
#' @param nulltype type of null distribution assumed in computing the probability of gene differential expression. 0 is the theoretical null N(0,1), 1 is maximum likelihood estimation.
#' @param df the degrees of freedom used for estimating marginal distrution.
#' @param pct the percentage of genes exculded from fiting the two-group mixture model.
#' @param plot plots desired.  0 gives no plots. 1 gives single plot showing the
#' histogram of zz and fitted densities f and p0*f0.
#' @param spar smoothing parameter used to fit a smoothing spline, typically
#' (but not necessarily) in (0,1].  The coefficient lambda of the integral
#' of the squared  second derivative in the fit (penalized log likelihood)
#' criterion is a monotone function of ‘spar’, see the details from \code{help(smooth.spline)}
#' @param estBaseline compute the dispersion parameter only using the baseline transcriptome
#' @param estSize if True, size parameter is estiamted from each bin; if False,
#' dispersion parameter is estiamted from each bin.
#' 
#' 
#' @return 'iDEG' produces a list containing the following elements:
#'
#' \describe{
#'   \item{results}{a table iDEG result for each gene. The first two columns are
#' the gene epxression values of the two transcriptomes provided by the user.
#' The thrid column is the local false discovery rate, which provides the
#' probability of a gene being differentially expresseed. The fourth column is
#' the statistic used to compute the local false discovery rate, and can be used
#' as an effect size.}
#' \item{sizeHat}{When the assumptioin of constant dispersion across genes is
#' made, this is an single estimate of the common dispersion. When the
#' assumptioin of non-constant is made, this is a vector of estimates for the
#' dispersion parameter of each gene.}
#' }
#' 
#' @examples
#' set.seed(1)
#' exp_mean1 <- rexp(20000, 1/500) + 1
#' exp_mean2 <- exp_mean1
#' exp_mean2[1:100] <- exp_mean2[1:100] * 10
#' transcriptome1 <- rnbinom(n = length(exp_mean1), size = 60, mu = exp_mean1)
#' transcriptome2 <- rnbinom(n = length(exp_mean2), size = 60, mu = exp_mean2)
#' res <- iDEG(transcriptome1,transcriptome2)
#'
#' @export
 
#####
iDEG <- function (baseline, case, normalization = F,
                  dataDistribution = c('NB','Poisson'), numBin = 100,
                  rankBaseline = T,
                  estBaseline = F,
                  estSize = F,
                  spar = NULL,
                  plot = 0,
                  constDisp = T,  df = 7, nulltype = 1, pct = .0001 ){
    
    ## verify that the two transcriptomes contain the same number of genes.
    message('make sure sample 1 is baseline')
    if( !is.vector(baseline) | !is.vector(case) ) stop('The two provied transcriptomes need to be numerical vectors, and the name of each element in the vectors should be the gene name corresponding to that numercial value')
    if( length(baseline) != length(case) ) stop('The two provided transcriptomes contain different numbers of genes.')
    if( any(names(baseline) != names(case)) ) stop('The genes of the two provided transcriptomes are not in the same order. Please check the gene names.')

    dataDistribution <- match.arg(dataDistribution)
    
    ## combine the two transcriptomes to make a data frame
    data1 <- data.frame(baseline,case, check.rows = T)
    
    ## record the total number of genes
    totNumGenes <- length(baseline)
    
    ## remove the genes with both values being zero
    ind0 <- rowSums(data1) == 0 # these genes will be claimed as non-DEG
    data_test <- data1[!ind0,]
    if(dataDistribution == 'NB'){ 
        iDEG.NB(dataTest= data_test, data1 = data1, ind0 = ind0,
                normalization = normalization,
                numBin = numBin,estSize =estSize,
                rankBaseline =rankBaseline, estBaseline = estBaseline,
                constDisp = constDisp, spar = spar, df = df,
                nulltype = nulltype, pct = pct, plot = plot)
    }else if(dataDistribution == 'Poisson'){
        iDEG.Pois(dataTest = data_test, data1 = data1,
                  ind0 = ind0, normalization = normalization,
                  df = df, nulltype = nulltype, pct = pct, plot = plot)
     }
}



