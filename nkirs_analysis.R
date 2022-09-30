## ----knitrSetup---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=TRUE, message=FALSE, warning=FALSE, fig.width=4.5, fig.height=4.5, cache=FALSE)
knitr::knit_hooks$set(GPs=function(before, options, envir){
  if (before){
    cex=1.5
    # graphical parameters
    fgcol = 'white'
    bgcol = 'black'
    par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
        cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
        mar=c(0.5,2.2,2,0), bty='l', oma=c(0,0,2,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
  }})
knitr::opts_chunk$set(echo = FALSE, fig.height = 4, fig.width = 4, GPs=TRUE, cache=FALSE, cache.lazy=FALSE, eval=TRUE)
path = Sys.getenv('PATH')
path = Sys.setenv('PATH'=paste(path, '/home/rstudio/.local/bin', sep=':'))
set.seed(555)


## ----dataSetup, eval=TRUE-----------------------------------------------------
# install the latest versions of the packages to perform these analyses.
devtools::install_github('simonvandekar/pbj', ref='inference')
#devtools::install_github('statimagcoll/NIsim')

### LIBRARIES ###
library(RNifti) # Nifti I/O
library(parallel) # mclapply
library(mmand) # spatial cluster functions
library(fslr) # imaging tools
library(progress) # progress bar
library(pbj) # pbj package
library(PDQutils) # edgeworth expansion stuff
# library(NIsim) # simulation tools
library(papayar) # image viewer
library(splines) # ns
library(magrittr) # %>%

# number of cores for parallel things
ncores = 16


### LOAD IN DATA FROM DROPBOX ###
dbimagedir = '/media/disk2/pbj/data/rockland/neuroimaging'
# maskfile = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
# new mask file (created May 12, 2021)
maskfile = "/media/disk2/pbj/data/rockland/neuroimaging/overlap_mask_2mm.nii.gz"
dbdatafile = '/media/disk2/pbj/data/rockland/demographic/RocklandBehavioral.csv'
datafile = '/media/disk2/pbj/pbj_ftest/nkirs_bootstrap_results.rdata'
templatefile = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz'

# # creates slab mask
# mask = readNifti(maskfile)
# slabs = round(dim(mask)[3]/2) + -3:3
# sum(mask[,,slabs])
# mask[,,-slabs] = 0
# writeNifti(mask, '/media/disk2/pbj/data/rockland/neuroimaging/MNI152_T1_2mm_brain_mask_slab.nii.gz')



# load in data and get directories
dat = read.csv(dbdatafile)
dat$dir = file.path(dbimagedir, dat$AnonymizedID, dat$IDandSession)
#dat$dir[which(!file.exists(dat$dir))]
dat = dat[ file.exists(dat$dir), ]
# some subjects have two image folders, this just grabs the first one.
dat$files = file.path(sapply(dat$dir, function(dir) list.files(dir, pattern='*', full.names=TRUE)[1]), 'GRAY_MNORM/mwp1t1.nii.gz')
dat$files2mm = gsub('.nii.gz', '_2mm.nii.gz', dat$files)
dat$files2mmsm4 = gsub('.nii.gz', '_sm8.nii.gz', dat$files2mm)


if(!all(file.exists(dat$files2mmsm4))){
  ## CREATE DOWNSAMPLED DATA
  invisible(mcmapply(flirt, infile=dat$files, outfile=dat$files2mm, opts = '-applyxfm', MoreArgs=list(reffile=templatefile, retimg=FALSE), mc.cores = ncores ))
  ## SMOOTH DOWNSAMPLED DATA
  # convert FWHM to sigma
  sigma = 8/2.355
  invisible(mcmapply(susan, file=dat$files2mm, outfile=dat$files2mmsm4, MoreArgs=list(sigma=sigma, dimg='3', n_usans='0'), mc.cores=ncores)) 
  # sigma = 8: using 8mm smoothing (?)
}

# view files
# papayar::papaya(c(template, dat$files2mm[1], dat$files2mmsm4[1]))
# set downsampled smoothed files as the main file
#papayar::papaya(c(maskfile, dat$files2mm[1], dat$files2mmsm4[1]))
dat$files = dat$files2mmsm4

# rename some variables and do some data curation
dat$age = dat$CalculatedAge
dat$sex = dat$WhatIsYourSex
dat$race = ifelse(!dat$WhatIsYourRace %in% c('White', 'Black'), 'Other', dat$WhatIsYourRace)
dat$pdf = file.path(sapply(dat$dir, function(dir) list.files(dir, pattern='*', full.names=TRUE)[1]), 'PDF/catreport_t1.pdf')

## CREATE MASK
# if(!file.exists(maskfile)){
#   imgs = readNifti(dat$files)
#   mask = imgs[[1]]
#   mask[,,] = 0
#   imgs = apply(simplify2array(imgs), 1:3, function(v) sum(v>0))
#   mask[sum(imgs)==nrow(dat)] = 1
#   writeNifti(mask, file=maskfile)
# }



## ---- model-------------------------------------------------------------------
# Testing the nonlinearity of age effect
# lmfull = paste0(" ~ sex + race + ns(age, df=4)" )
# lmred = paste0(" ~ sex + race + age" )

# Testing the interaction between nonlinear age x sex
lmfull = paste0(" ~ race + ns(age, df=4) * sex" )
lmred = paste0(" ~ sex + race + ns(age, df = 4)" )
dat = dat[apply(!is.na(dat[,all.vars(as.formula(lmfull))]), 1, all), ]


## -----------------------------------------------------------------------------

# image(paramStatMap)
# lmPBJ -- fits the model and computes the statistical image for the covariate of interest
paramStatMap = lmPBJ(dat$files, form=lmfull, formred=lmred, mask=maskfile, template = templatefile, data=dat, transform = 'none', HC3 = TRUE, robust = FALSE)

# Visualize the chi square image
# statmapFile = paste0(tempfile(), '.nii.gz') # generate a random file name
# writeNifti(stat.statMap(paramStatMap), file=statmapFile) # write the statistical image
# papayar::papaya(c(template, statmapFile )) # view the statistical image



## -----------------------------------------------------------------------------
# lmPBJ -- fits the model and computes the statistical image for the covariate of interest
## using robust test statistics (`robust = TRUE`)
robustStatMap <- lmPBJ(dat$files, form=lmfull, formred=lmred, mask=maskfile, template = templatefile, data=dat, transform = 'none', HC3 = TRUE, robust = TRUE)


# Visualize the chi square image
# statmapFile = paste0(tempfile(), '.nii.gz') # generate a random file name
# writeNifti(stat.statMap(robustStatMap), file=statmapFile) # write the statistical image
# papayar::papaya(c(template, statmapFile )) # view the statistical image


## ---- pbjInferenceSetup-------------------------------------------------------
# cluster extent inference
## cluster forming threshold
cft = qchisq(0.01, df = paramStatMap$sqrtSigma$df, lower.tail = FALSE)
mask = RNifti::readNifti(maskfile)
# number of bootstraps (or permutations) to run
nboot=1000

# This function is used to compute local maxima and cluster extents within each bootstrap/permutation
statistic = function(stat, rois=FALSE, mask, thr){
  c(list(maxima(stat, rois=rois)), cluster(stat, mask=mask, thr=thr, rois=rois) )
}


## -----------------------------------------------------------------------------
  ## wild bootstraps with rademacher
  param.wild.time = system.time(invisible(capture.output(pbjinf.param.wild <- pbjInference(paramStatMap, statistic = statistic, nboot = nboot, runMode = 'cdf', method='t', mask = paramStatMap$mask, thr = cft))))

## permutations
  param.permu.time = system.time(invisible(capture.output(pbjinf.param.permu <- pbjInference(paramStatMap, statistic = statistic, nboot = nboot, runMode = 'cdf', method='permutation', mask = paramStatMap$mask, thr = cft) )))


## -----------------------------------------------------------------------------
## wild bootstraps with rademacher
robust.wild.time = system.time(invisible(capture.output(pbjinf.robust.wild <- pbjInference(robustStatMap, statistic = statistic, nboot = nboot, runMode = 'cdf', method='t', mask = robustStatMap$mask, thr = cft))))

## permutations
robust.permu.time = system.time(invisible(capture.output(pbjinf.robust.permu <- pbjInference(robustStatMap, statistic = statistic, nboot = nboot, runMode = 'cdf', method='permutation', mask = robustStatMap$mask, thr = cft) )))

save.image(datafile)


## ----paramClusterTable, comment=NA, eval=TRUE---------------------------------
load(datafile)
Param.Table = data.frame('Cluster Extent' = c(pbjinf.param.wild$obsStat[[2]]), 
                         'Unadjusted p-value - boot' = (1-pbjinf.param.wild$margCDF[[2]](c(pbjinf.param.wild$obsStat[[2]]))), 
                         'Unadjusted p-value - perm' = (1-pbjinf.param.permu$margCDF[[2]](c(pbjinf.param.permu$obsStat[[2]]))),
                         'FWER p-value - boot' = (1-pbjinf.param.wild$globCDF[[2]](c(pbjinf.param.wild$obsStat[[2]]))),
                         'FWER p-value - perm' = (1-pbjinf.param.permu$globCDF[[2]](c(pbjinf.param.permu$obsStat[[2]]))),
                         check.names=FALSE )

knitr::kable(Param.Table[order(Param.Table$`Cluster Extent`, decreasing = TRUE),][1:10,], row.names = FALSE,
             digits=3, booktabs=T)
print(knitr::kable(Param.Table[order(Param.Table$`Cluster Extent`, decreasing = TRUE),][1:10,], format = 'latex', row.names = FALSE,
             digits=3, booktabs=T))


## ----robustClusterTable, comment=NA, eval=TRUE--------------------------------
Robust.Table = data.frame('Cluster Extent' = c(pbjinf.robust.wild$obsStat[[2]]), 
                         'Unadjusted p-value - boot' = (1-pbjinf.robust.wild$margCDF[[2]](c(pbjinf.robust.wild$obsStat[[2]]))), 
                         'Unadjusted p-value - perm' = (1-pbjinf.robust.permu$margCDF[[2]](c(pbjinf.robust.permu$obsStat[[2]]))),
                         'FWER p-value - boot' = (1-pbjinf.robust.wild$globCDF[[2]](c(pbjinf.robust.wild$obsStat[[2]]))), 
                         'FWER p-value - perm' = (1-pbjinf.robust.permu$globCDF[[2]](c(pbjinf.robust.permu$obsStat[[2]]))),
                         check.names=FALSE )

knitr::kable(Robust.Table[order(Robust.Table$`Cluster Extent`, decreasing = TRUE),][1:10,], format = 'html', row.names = FALSE,
             digits=3, booktabs=T)
print(knitr::kable(Robust.Table[order(Robust.Table$`Cluster Extent`, decreasing = TRUE),][1:10,], format = 'latex', row.names = FALSE,
             digits=3, booktabs=T))


## ----figuresetup, eval=TRUE---------------------------------------------------
# read in template
template = readNifti(templatefile)
# PARAMETIC STATISTICS
displayParamBoot <- displayParamPerm <- stat.statMap(paramStatMap)

# 1. (wild) bootstrap parametric statistic using cluster extent inference
  # get the FWER-adjusted p-values
  pvals = round(-log10( 1 - pbjinf.param.wild$globCDF[[2]](c(pbjinf.param.wild$obsStat[[2]]) ) + 0.0001), 3)
  # getting the ROI image, which indexes what cluster each voxel belongs to
  pvalimg = pbjinf.param.wild$ROIs[[2]]
  # This uses ROI image to index the p-value vector and replace voxel values in the image with their adjusted p-value
  pvalimg[ pvalimg>0] = pvals[pvalimg]
  # This thresholds the unadjusted chi-square statistical image by the adjusted p-value for the cluster
  displayParamBoot[ pvalimg< (-log10(0.05)) ] = 0
  
# 2. permutation parametric statistic
  # get the FWER-adjusted p-values (for each cluster)
  pvals = -log10( 1 - pbjinf.param.permu$globCDF[[2]](c(pbjinf.param.permu$obsStat[[2]]) ) + 0.0001) %>% round(3)
  # getting the ROI image, which indexes what cluster each voxel belongs to
  pvalimg = pbjinf.param.permu$ROIs[[2]]
  # This uses ROI image to index the p-value vector and replace voxel values in the image with their adjusted p-value
  pvalimg[ pvalimg>0] = pvals[pvalimg]
  # bootstrap parametric statistic using cluster extent inference
  # This thresholds the unadjusted chi-square statistical image by the adjusted p-value for the cluster
  displayParamPerm[ pvalimg< (-log10(0.05)) ] = 0

# ROBUST STATISTICS
displayRobustBoot <- displayRobustPerm <- stat.statMap(robustStatMap) 

# 1. (wild) bootstrap robust statistic using cluster extent inference
  # get the FWER-adjusted p-values
  pvals = -log10( 1 - pbjinf.robust.wild$globCDF[[2]](c(pbjinf.robust.wild$obsStat[[2]]) ) + 0.0001) %>% round(3)
  
  # getting the ROI image, which indexes what cluster each voxel belongs to
  pvalimg = pbjinf.robust.wild$ROIs[[2]]
  
  # This uses ROI image to index the p-value vector and replace voxel values in the image with their adjusted p-value
  pvalimg[ pvalimg>0] = pvals[pvalimg]
  
  # This thresholds the unadjusted chi-square statistical image by the adjusted p-value for the cluster
  displayRobustBoot[ pvalimg< (-log10(0.05)) ] = 0
  
# 2. permutation robust statistic (using cluster extent inference)
  # get the FWER-adjusted p-values (for each cluster)
  pvals = -log10( 1 - pbjinf.robust.permu$globCDF[[2]](c(pbjinf.robust.permu$obsStat[[2]]) ) + 0.0001) %>% round(3)
  # getting the ROI image, which indexes what cluster each voxel belongs to
  pvalimg = pbjinf.robust.permu$ROIs[[2]]
  # This uses ROI image to index the p-value vector and replace voxel values in the image with their adjusted p-value
  pvalimg[ pvalimg>0] = pvals[pvalimg]
  # bootstrap parametric statistic using cluster extent inference
  # This thresholds the unadjusted chi-square statistical image by the adjusted p-value for the cluster
  displayRobustPerm[ pvalimg< (-log10(0.05)) ] = 0
  
  
  #displayRobustPerm[,,] = displayRobustPerm[fullRows,fullCols,]
  #displayRobustBoot[,,] = displayRobustBoot[fullRows,fullCols,]
  #displayParamBoot$ = displayParamBoot[fullRows,fullCols,]


## ---- height=4.5, width=2.5, eval=TRUE----------------------------------------
# color bar function
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, round(ticks, 2), las=1, cex.axis=cex*0.7, font=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# choose upper threshold that looks good
threshs = c(-log10(0.01), 26)

# display figure "figures for evidence"
# slices that have the mask in them
# which(!apply(mask==0, 3, all))
# a manually selected subset of that
slices = 10:60

# create this directory
outputdir = '/media/disk2/pbj/pbj_ftest/nkirs_images'

for (slice in slices){
  fname = file.path(outputdir, 'images', paste0('slice', slice, '.png') )
  dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
  png(filename = fname, height=4, width=4, units = 'in', res = 300)
  cex=1.5
  # graphical parameters
  fgcol = 'white'
  bgcol = 'black'
  par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
      cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
      mar=c(0,2.2,2.2,0), bty='l', oma=c(0,0,0,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
      layout(cbind(matrix(1:4, nrow=2, byrow=TRUE) %x% matrix(1, nrow=3, ncol=3) , 5))
      # display parametric bootstrap statistic
      image(displayParamBoot, template, thresh=(-log10(0.01)), index=slice, cex=cex)
      mtext('Parametric', side=2, cex = 0.8*cex, font = 2)
      mtext('Bootstrap', side = 3, cex = 0.8*cex, font = 2)
      # display parametric permutation statistic
      image(displayParamPerm, template, thresh=(-log10(0.01)), index=slice, cex=cex) 
      mtext('Permutation', side = 3, cex = 0.8*cex, font=2)
      # display Robust bootstrap statistic
      image(displayRobustBoot, template, thresh=(-log10(0.01)), index=slice, cex=cex) 
      mtext('Robust', side=2, cex = 0.8*cex, font = 2 )
      # display robust permutation statistic
      image(displayRobustPerm, template, thresh=(-log10(0.01)), index=slice, cex=cex) 
      
      # main title
      #mtext('Probability', side=3, outer = TRUE, cex=1*cex, font=2)
        
      par(mar=c(10,2,10,0.5))
      color.bar(pbj:::redyellow(64), min=threshs[1], max=threshs[2], nticks=4)
      dev.off()
}

