#' Generates data for simulation
#'
#' This function loads images specified in filelist and adds signal to the
#'  images. The amount and structure of signal is determined by betaimg and X.
#'  The output images are equal the images specified by filelist plus betaimg
#'  times the column of X that is not in Xred. This column is first
#'  residualized to Xred.
#' @param files a vector of .nii or .nii.gz images.
#' @param n4sim total number of measurements that we're simulating
#' @param nMeas Approximate average number of measurements per subject.
#' @param mask image mask where data exist.
#' @param method Method to generate data, either "synthetic" (i.e. multivariate normal) or "bootstrap."
#' @param lambda Simulated data are a convex coombination of normally distributed data with bootstrapped data. lambda=1 is fully bootstrapped data, lambda=0 is fully synthetic.
#' @param outfiles a vector of images to save the output.
#' @keywords  null power simulation
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom RNifti writeNifti
#' @export
# @examples
genSimData = function(files, outfiles=NULL, n4sim, nMeas=NULL, mask=NULL, method=c('bootstrap', 'synthetic'), lambda=0.5 ){

    if(tolower(method[1])=='synthetic'){
      # THIS CODE SIMULATES DATA FOR THE FULL SAMPLE, if n in simulation is larger than n in sample then it will create problems
      # contains residuals for entire study
      cov = if(length(files)==1) readRDS(files[1]) else files
      if(!is.list(cov)){
        y = (matrix(rnorm(nrow(cov)*ncol(cov)), nrow=ncol(cov), ncol=nrow(cov))  %*% cov)/sqrt(nrow(cov))
        rm(cov)
        temp = if(is.character(mask)) readNifti(mask) else mask
        trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
        RNifti::writeNifti(temp, outfiles[ind])
        })
        rm(y)
      } else {
        # This is a hybrid of a bootstrap and parametric simulation
        # generate y using sample residuals, coefficients, and random effects
        n = nrow(cov$resids)
        nRE = nrow(cov$ranefs)
        N4sim = round(n4sim * nMeas)
        id = rep(1:n4sim, ceiling(N4sim/n4sim))[1:N4sim]
        samp = sample(n, N4sim, replace=TRUE)
        reSamp = sample(nRE, n4sim, replace=TRUE)
        # currently written only for random intercept
        REs = matrix(rnorm(n4sim*nRE)/sqrt(nRE), nrow=n4sim) %*% cov$ranefs[,1,]
        epsilon = matrix(rnorm(N4sim * n)/sqrt(n), nrow=N4sim) %*% cov$resids
        # convex combination of bootstrap resample and parametric resample
        # also assumes random intercept model only
        REs = (1-lambda)* REs + lambda * cov$ranefs[reSamp,1,]
        y = REs[id,] + (1-lambda) * epsilon + lambda * cov$resids[samp,]

        # write out Nifti images
        temp = if(is.character(mask)) readNifti(mask) else mask
        trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
        RNifti::writeNifti(temp, outfiles[ind])
        })
        # return simulated data if requested
        list(outfiles=outfiles, simdata=y, id=id)
        }
    } else if(tolower(method[1])=='bootstrap'){
      # the bootstrapping has already been performed in the simulation Setup phase
      result = file.copy(files, outfiles)
    } else{
      stop('genSimData method is not correctly specified.')
  }
}
