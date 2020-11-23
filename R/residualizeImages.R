#' Prepares Simulation Data for Bootstrapping
#'
#' Residualizes the images in files to the model form and writes the output to
#' outfiles optionally smoothes with sm (in mm FWHM) using susan prior to
#' residualizing if smoutfiles is specified smoothed images are saved into that
#' directory. Used to run simulations to assess power and type 1 error for
#' papers. This creates images that are residualized to the covariates which
#' can then be bootstrapped to generate a sample where there is the potential
#' for heteroskedasticity/nonexchangeability, but where the covariates are
#' unassociated with the mean of the outcome.
#' @param files Character vector of subject images to be modeled as an outcome
#'  variable.
#' @param form mgcv or lm style formula.
#' @param dat Data frame containing covariates used by form.
#' @param mask Character giving location of mask image.
#' @param outfiles Character vector of residual output images.
#' @param outrds Character vector of RDS file to save residuals as a matrix.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @return No returned value. This functions saves out nifti images files after residualizing to the model specified by form and dat. The residuals of files are saved as the corresponding element in outfiles.
#' @keywords power simulation, parametric bootstrap, type 1 error simulations, null simulations
#' @importFrom RNifti writeNifti
#' @importFrom parallel mclapply
#' @importFrom stats model.matrix
#' @export
# @examples
residualizeImages = function(files, form, dat, mask, outfiles=NULL, outrds=NULL, mc.cores=getOption("mc.cores", 2L)){

  cat('loading images.\n')
  y = simplify2array(mclapply(files, readNifti, mc.cores=mc.cores))

  # run linear model to get residuals
  if(is.character(mask)){
    mask = RNifti::readNifti(mask)
  }
  y = t(apply(y, length(dim(y)), function(x) x[mask==1]))
  X = model.matrix(form, data=dat)
  cat('regressing out covariates.\n')
  y = qr.resid(qr(X), y)
  invisible(sapply(dirname(outfiles), dir.create, showWarnings=FALSE) )

  if(!is.null(outrds)){
    saveRDS(y, outfiles)
  }
  if(!is.null(outfiles) ) {
    # save out images
    temp = mask
    invisible(lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
    RNifti::writeNifti(temp, outfiles[ind])
    }))
  }
  outfiles
}
