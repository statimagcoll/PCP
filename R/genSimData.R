#' Generates data for simulation
#'
#' This function loads images specified in filelist and adds signal to the
#'  images. The amount and structure of signal is determined by betaimg and X.
#'  The output images are equal the images specified by filelist plus betaimg
#'  times the column of X that is not in Xred. This column is first
#'  residualized to Xred.
#' @param files a vector of .nii or .nii.gz images.
#' @param mask image mask where data exist.
#' @param method Method to generate data, either "synthetic" (i.e. multivariate normal) or "bootstrap."
#' @param betaimg a parameter image that describes the association between X
#'  and the desired relationship images in filelist. The units of this are standardized so that it
#'  is like effect size. Currently, not supported.
#' @param outfiles a vector of images to save the output.
#' @keywords  null power simulation
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom RNifti writeNifti
#' @export
# @examples
genSimData = function(files, outfiles=NULL, betaimg=NULL, mask=NULL, method=c('bootstrap', 'synthetic') ){

  if(is.null(betaimg)) {
    if(tolower(method[1])=='synthetic'){
      # contains residuals for entire study
    cov = readRDS(file.path(dirname(files[1]), 'residuals.RDS'))
    y = (matrix(rnorm(nrow(cov)*length(files)), nrow=length(files), ncol=nrow(cov))  %*% cov)/sqrt(nrow(cov))
    rm(cov)
    temp = if(is.character(mask)) readNifti(mask) else mask
    trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
    RNifti::writeNifti(temp, outfiles[ind])
    })
    rm(y)
    } else if(tolower(method[1])=='bootstrap'){
      # the bootstrapping has already been performed in the simulation Setup phase
      result = file.copy(files, outfiles)
    } else{
      stop('genSimData method is not correctly specified.')
    }
  }
}
