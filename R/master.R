#' Sets up a Neuroimaging bootstrap-based simulation given a set of images and covariates.
#'
#' The images can be images from a real data set.
#' The bootstrap-based simulation conditions on the distribution of the sample by drawing subsets with replacement from the sample.
#' @param images String vector containing paths to the images used for the simulation.
#' @param data A dataframe with number of rows equal to the length of the images variable with corresponding covariates.
#' @param outdir A directory to save the output files that are used for the simulations.
#' @param nsim Number of simulations to setup for each sample size.
#' @param ns Sample sizes to evaluate for each simulation. Should be less than the number of images.
#' @param mask If performing simulations under an alternative, signal can be added to the images within the mask.
#' @param rs vector of radii for signal spheres.
#' @param betas vector of parameters for signal spheres. Signal is constant throughout the sphere.
#' @return Returns a matrix of directories where simulation setup files are stored.
#' @importFrom RNifti readNifti writeNifti
#' @importFrom foreach foreach %dopar%
#' @export
# prepare the output directories
# All the randomization happens within this loop
simSetup = function(images, data, outdir, nsim=1000, ns=c(50,100, 200, 400), mask=NULL, rs=8, betas=rep(0, length(rs)) ){
  sims = expand.grid(sim=1:nsim, n=ns, simdir=NA)
  data$images = images
  if(any(betas!=0) & is.null(mask)) stop('mask must be provided if simulating under an alternative.')
  foreach(simind=1:nrow(sims), .combine=list) %dopar% {
    simdir = sims$simdir[simind] = file.path(outdir, paste0('sim', sims[simind,'sim']), paste0('n', n) )
    n = sims[simind, 'n']
    dir.create(file.path(simdir, paste0('n',n)), showWarnings=FALSE, recursive = TRUE)
    unlink(file.path(simdir, '*.nii.gz'), recursive=TRUE)
    # create a sphere where there is signal within the gray matter
    # mask it with the (gray matter) mask.
    if(any(betas!=0)) parameterImage(mask, parameterImage = file.path(outdir, 'signal.nii.gz'), rs, betas)
    # create random sample from demographics and roi data
    tempdata = data[sample.int(nrow(data), n, replace=TRUE), ]
    saveRDS(tempdata, file=file.path(simdir,'data.rds' ) )
  }# end for(sim)
  sims
}


  #' Creates parameter image in random locations for simulations
  #'
  #' Uses the mask variable to create length(rs) spheres of size rs with values equal to betas within the mask.
  #' A random voxel is chosen within the mask and a sphere is placed at that location.
  #' The spheres are masked by the mask image so that no parameter values exist outside the mask.
  #' @param mask String or niftiImage object indication the area to select sphere locations from.
  #' @param parameterImage output nifti for the parameter image.
  #' @param rs vector of radii for signal spheres.
  #' @param betas vector of parameters for signal spheres. Signal is constant throughout the sphere.
  #' Model parameter values are set to params = betas * sd(y)/sd(x).
  #' @return Returns the parameter image after writing it to file.path(outdir, 'signal.nii.gz').
  #' @importFrom RNifti readNifti writeNifti
  #' @export
  parameterImage = function(mask, parameterImage, rs, betas){
    if(is.character(mask)) mask = readNifti(mask)
    outfile = mask
    inds = which(mask==1, arr.ind=TRUE)
    # random location in gray matter mask for center voxel
    centers = inds[sample.int(nrow(inds), length(rs)),, drop=FALSE]
    outfile[,,] = 0
    for(rind in 1:length(rs)){
      r = rs[rind]
      center = centers[rind,]
      inds = as.matrix(expand.grid(seq(-r,r), seq(-r,r), seq(-r,r)))
      inds = inds[ sqrt(rowSums(inds^2))<=r,]
      inds = sweep(inds, 2, center, '+') # sphere around center voxel
      # make sure all voxels are in the image
      inds = inds[apply(sapply(1:3, function(x) inds[,x]>0 & inds[,x]<=dim(outfile)[x]), 1, all), ]
      # It is possible that a small cluster can be captured inside a big cluster and then it won't exist in that simulation
      outfile[ inds ] = betas[rind]
      # signal only in the mask
    } # end for(rind)
    outfile = outfile * mask
    writeNifti(outfile, file=parameterImage)
    outfile
  }

  #' Interface to command line tool for dropbox
  #'
  #' @param cmd String indicating command to run.
  #' @param ... arguments passed to command.
  #' @param ncores number of vectors of commands.
  #' @importFrom parallel mclapply
  #' @export
  dbxcli = function(cmd, ..., ncores=1){
    result = parallel::mclapply(paste('dbxcli', cmd, do.call(paste, list(...) ) ), system, mc.cores=ncores )
  }



  #' Configure AWS
  #'
  #' Configures AWS using access key ID and secret access key provided in csv by user.
  #' Need to run on master. This sets the default access key and id.
  #' to get access key, click on name in top right of the management page, click security credentials, click access keys, click create access key.
  #' @param keycsv character path to csv file obtained from AWS containing the id and key.
  #' @param profile character username to setup profile.
  #' @param region character amazon region to use.
  #' @param output character amazon default output from aws command line interface tool.
  #' @export
  #' @importFrom utils read.table
  configureAWS = function(keycsv, profile='default', region='us-east-2', output='json'){
    dir.create('/home/rstudio/.aws', showWarnings = FALSE)
    if(!is.null(keycsv)){
      suppressWarnings(key <- read.table(keycsv, stringsAsFactors = FALSE, sep=',', header = FALSE))
      id = gsub('.*=', '', key[1,1])
      key = gsub('.*=', '', key[2,1])

      fileConn<-file("~/.aws/credentials")
      writeLines(c(profile, id , key), fileConn)
      close(fileConn)

      fileConn<-file("~/.aws/config")
      writeLines(c(profile, region, output), fileConn)
      close(fileConn)
    }
  }

  #' List memory usage of all objects
  #'
  #' @param units character what units to use, passed to format for object.size.
  #' @export
  #' @importFrom utils object.size
  memoryUse = function(units='MiB'){
    sort( sapply(ls(),function(x){format(object.size(get(x)), units=units)}))
  }
