## ----setup, include=FALSE-----------------------------------------------------
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
knitr::opts_chunk$set(echo = FALSE, fig.height = 4, fig.width = 4, GPs=TRUE, cache=TRUE, cache.lazy=FALSE)
path = Sys.getenv('PATH')
path = Sys.setenv('PATH'=paste(path, '/home/rstudio/.local/bin', sep=':'))
set.seed(555)


## ----simsetup-----------------------------------------------------------------
# install the latest versions of the packages to perform these analyses.
devtools::install_github('simonvandekar/pbj', ref='ftest')
devtools::install_github('statimagcoll/NIsim', ref='master')

### LIBRARIES ###
library(RNifti)
library(parallel)
library(splines)
library(progress)
library(pbj)
library(NIsim)

### DATA LOCATION ###
datadir = '/media/disk2/pbj/data'
# will be created by downloadABIDE
dbimagedir = file.path(datadir, 'abide/neuroimaging/cpac/alff')
templatefile = file.path(datadir, 'abide/neuroimaging/MNI152_T1_3mm.nii.gz')
# created by simSetup... Not used in this analysis
dbresimagedir = file.path(datadir, 'abide/neuroimaging/cpac/alff_res')
# created below
maskfile = file.path(datadir, 'abide/neuroimaging/cpac/mask.nii.gz')
#dbdatafile = '~/pbj/data/abide/demographic/n1035_phenotypic_20190509.rds'
#dbimagedir = '~/pbj/data/abide/neuroimaging/cpac/alff_cropped/'
#dbresimagedir = '~/pbj/data/abide/neuroimaging/cpac/alff_cropped_res/'
#maskfile = '~/pbj/data/abide/neuroimaging/cpac/cropped_n1035_mask.nii.gz'

# load in data and get directories
dat = downloadABIDE(datadir)
# some data curation is done on download. Subset all pass on visual inspection and making sex and dx factors, remove rows with no file name
dat$imgname = paste(dat$file_id, 'alff.nii.gz', sep='_')
dat$files = file.path(dbimagedir, dat$imgname)


  imgs = simplify2array(readNifti(dat$files) )
  # choose
  # number of people with zeros at that location
  inds=numeric()
  ids = c()
  subids = dat$sub_id
  # number of voxels with no zeros
  nnzero = 0
  # iteratively removes subjects who will increase the mask the largest
  while(nnzero<30000){
    voxSums = rowSums(imgs==0, dims=3)
    tab = as.data.frame(table(voxSums))
    nnzero=tab[1,2]
    # number of unique voxels for each subject
    uniquevox = apply(imgs, 4, function(img) sum(img==0 & voxSums==1) )
    # number of subjects to remove based on those subjects decreasing the amount of unique zero voxels by 50%
    #ind = which(cumsum(sort(uniquevox[ uniquevox>0], decreasing = TRUE))/tab$Freq[2]>0.5)[1]
    #inds = which(uniquevox>=sort(uniquevox, decreasing = TRUE)[ind])
    inds = which.max(uniquevox)
    cat('\nIteration:\nsubject removed: ', paste(subids[inds], collapse=', '), '\nmask size is now ', nnzero+sum(uniquevox[inds]), '\nNumber of voxels added:', sum(uniquevox[inds]) )
    imgs = imgs[,,,-inds]
    subids = subids[-inds]
    nnzero=nnzero+sum(uniquevox[inds])
  }
  nexcluded = nrow(dat) - length(subids)
  cat('\nExcluded', nexcluded, 'subjects.')
  #tab$cumvox = rev(cumsum(rev(tab$Freq)))
  dat = dat[ dat$sub_id %in% subids,]
  
  # now create the mask
  mask = readNifti(dat$files[1])
  imgs = simplify2array(readNifti(dat$files) )
  # mask actually still has pretty low coverage
  mask[,,] = 0
  mask[,,] = apply(imgs>0, 1:3, all)
  writeNifti(mask, maskfile)
  nvox = sum(mask)
  rm(imgs)


## ---- simconfig---------------------------------------------------------------
### SIMULATION PARAMETERS ###
simConfig = list(
  # vector of sample sizes to simulate
  ns = 25 * 2^(0:5),
  # number of simulations to run
  nsim=10000,
  nboot=200,
  # cluster forming thresholds
  cfts.s = c(0.1, 0.2, 0.25, 0.4),
  cfts.p = c(0.05, 0.01, 0.001),
  
  # radius for spheres of signal.
  rs=c(8),
  
  #### MODEL FORMULAS FOR SIMULATIONS ####
  formres = paste0(" ~ dx_group + sex + ns(func_mean_fd, df=10) + ns(age_at_scan, df=10)" ),
  # need age_at_scan in both models for testing nonlinear functions
  form = paste0(" ~ sex + func_mean_fd + dx_group + age_at_scan" ),
  formred = paste0(" ~ sex + func_mean_fd + dx_group"),
  #  weights for each subject. Can be a character vector
  W = c("func_mean_fd"),
  # where to put residuals
  resdir = dbresimagedir,
  # where to output results
  simdir = '~/temp',
  dat = dat,
  mask = maskfile,
  ncores = 24,
  method='bootstrap'
)
simConfig$betas = rep(0, length(simConfig$rs))
simConfig$output = file.path(datadir, '../pbj_ftest', paste0('EST_reproducible_nsim', simConfig$nsim, '.rdata'))


## ----simulationFunctions------------------------------------------------------
simFunc = function(lmfull, lmred, mask, data, cfts.s, cfts.p, nboot, sim, seiIts=10){
  HC3RobustStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 'none', HC3 = TRUE )
  # t transform, classical, estimate covariance
  tStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 't', robust=FALSE, HC3=TRUE)
  
  # These probably don't work
  if(!sim %% seiIts){
    invisible(capture.output(seiResultsT <- pbjSEI(tStatmap, cfts.s=cfts.s, nboot = nboot, method='t')))
    invisible(capture.output(seiResultsRobust <- pbjSEI(HC3RobustStatmap, cfts.s=cfts.s, nboot = nboot, rboot=function(n){ (2*rbinom(n, size=1, prob=0.5)-1)}, method='t') ))
    invisible(capture.output(seiResultsTp <- pbjSEI(tStatmap, cfts.p=cfts.p, nboot = nboot, method='t')))
    invisible(capture.output(seiResultsRobustp <- pbjSEI(HC3RobustStatmap, cfts.p=cfts.p, nboot = nboot, rboot=function(n){ (2*rbinom(n, size=1, prob=0.5)-1)}, method='t')))
  } else {
    seiResultsT = NA
    seiResultsRobust = NA
    seiResultsTp = NA
    seiResultsRobustp = NA
  }
  
  out = list('tStatmap' = tStatmap$stat, 'HC3RobustStatmap'=HC3RobustStatmap$stat, 'tSEI'=seiResultsT, 'HC3RobustSEI'=seiResultsRobust, 'tSEIp'=seiResultsTp, 'HC3RobustSEIp'=seiResultsRobustp)
  gc()
  return(out)
}


#simdirs = simSetup(simConfig$dat$files, data=simConfig$dat, outdir=simConfig$simdir, nsim=simConfig$nsim, ns=simConfig$ns, mask=simConfig$mask, rs=simConfig$rs, betas=simConfig$betas )
#simtime = system.time(test <- simFunc(simConfig$form, simConfig$formred, simConfig$mask, readRDS(file.path(simdirs$simdir[200], 'data.rds')), 2, cfts.s = simConfig$cfts.s, nboot=simConfig$nboot, sim=simdirs$sim[200]))


## ----runSims, eval=TRUE-------------------------------------------------------
if(!file.exists(simConfig$output)){
  ### SETUP THE SIMULATION ANALYSIS ###
  # subsets dataset to all people who have the variables
  simConfig$dat = simConfig$dat[apply(!is.na(simConfig$dat[ ,c(all.vars(as.formula(simConfig$formres)), simConfig$W)]), 1, all), ]
  
  # setup the simulation output directories
  simdirs = simSetup(simConfig$dat$files, data=simConfig$dat, outdir=simConfig$simdir, nsim=simConfig$nsim, ns=simConfig$ns, mask=simConfig$mask, rs=simConfig$rs, betas=simConfig$betas )
  
  # test one simulation
  #time = system.time(test <- simFunc(simConfig$form, simConfig$formred, simConfig$mask, readRDS(file.path(simdirs$simdir[10], 'data.rds')), simConfig$nboot, simConfig$cfts.s) )
  
  # mix this up so that large sample simulations aren't all dropped on one "thread".
  simdirs = simdirs[sample(1:nrow(simdirs)),]
  
  results = runSim(simdirs=simdirs$simdir, sim=simdirs$sim, method=simConfig$method,
                   simfunc = simFunc, mask = simConfig$mask,
                   simfuncArgs = list(
                     lmfull= simConfig$form,
                     lmred = simConfig$formred,
                     mask = simConfig$mask,
                     cfts.s=simConfig$cfts.s,
                     cfts.p=simConfig$cfts.p,
                     nboot=simConfig$nboot), ncores = simConfig$ncores)
  
  dir.create(dirname(simConfig$output), showWarnings = FALSE, recursive = TRUE)
  # clean up files
  save.image(file=simConfig$output)
  unlink(list.files(tempdir(), full.names = TRUE))
  gc()
  stop('not an error. Finished simulations.')
}


## ---- figuresetup-------------------------------------------------------------
# also sets up the data frame with the output

# initialize output for loop
blankimg = mask = readNifti(maskfile)
blankimg[,,] = 0
template = readNifti(templatefile)

npns = c(length(simConfig$cfts.p), length(simConfig$cfts.s))
allimgout = data.frame(method=rep(c('p', 'S'), npns * length(simConfig$ns)), n=rep(simConfig$ns, sum(npns)), value = c(rep(simConfig$cfts.p, each=length(simConfig$ns)), rep(simConfig$cfts.s, each=length(simConfig$ns))), power=NA, mean=NA)
load(simConfig$output)
# fixes bad output from runSim... bug fixed in code, need to push to github.
#results = lapply(1:nrow(simdirs), function(ind) results[(1 +(ind-1)* 6):(1 +(ind-1)* 6 + 5)] )
#results = lapply(results, function(x){ names(x) = c('tStatmap', 'HC3RobustStatmap', 'tSEI', 'HC3RobustSEI', 'tSEIp', 'HC3RobustSEIp'); x})
simdirs$results = results
for(rowInd in 1:nrow(allimgout)){
  cat(rowInd, '\n')
  n = allimgout[rowInd, 'n']
  method = allimgout[rowInd, 'method']
  value = allimgout[rowInd, 'value']
  
  design = pbj::getDesign(simConfig$form, simConfig$formred, data=simConfig$dat)
  if(method=='p'){
    chisqValue = qchisq(value, df=design$df, lower.tail=FALSE)
  } else {
    #chisqValue = value^2*(n-ncol(design$X)) + design$df
    chisqValue = value^2*(n) + design$df
  }
  
  # Tstatmap
  blankimg[ mask==1] = rowMeans(simplify2array(lapply(simdirs$results[simdirs$n==n], function(x) x$tStatmap ))>=chisqValue)
  allimgout$tStatmapPower[rowInd] = list(blankimg)
  blankimg[,,] = 0
  temp = simplify2array(lapply(simdirs$results[simdirs$n==n], function(x) x$tStatmap ))
  temp[ is.infinite(temp)] = max(temp[is.finite(temp)])
  blankimg[ mask==1] = rowMeans(temp)
  allimgout$tStatmapMean[rowInd] = list(blankimg)
  rm(temp)
  #image(blankimg, template, thresh=chisqValue, index=slice, cex=cex*0.7)
  #title(paste(method, '=', value, ', n', '=', n, ', T-stat'))
  
  # Robust statmap
  blankimg[,,] = 0
  blankimg[ mask==1] = rowMeans(simplify2array(lapply(simdirs$results[simdirs$n==n], function(x) x$HC3RobustStatmap ))>=chisqValue)
  allimgout$robustStatmapPower[rowInd] = list(blankimg)
  blankimg[,,] = 0
  temp = simplify2array(lapply(simdirs$results[simdirs$n==n], function(x) x$HC3RobustStatmap ))
  temp[ is.infinite(temp)] = max(temp[is.finite(temp)])
  blankimg[ mask==1] = rowMeans(temp)
  allimgout$robustStatmapMean[rowInd] = list(blankimg)
  rm(temp)
  
  if(method=='S'){
    cftname = paste0('cft.s', value)
    # t type 1 error/power maps
    blankimg[,,] = 0 
    er = lapply(simdirs$results[simdirs$n==n], function(x){if(!is.na(x$tSEI[1])){ x$tSEI[[cftname]]$pmap> -log10(0.05)} else NA })
    er = er[ !is.na(er)]
    blankimg[,,] = rowMeans(simplify2array(er), dims=3)
    allimgout$tStatMapSEIpower[rowInd] = list(blankimg)
    # robust type 1 error/power maps
    er = lapply(simdirs$results[simdirs$n==n], function(x){if(!is.na(x$HC3RobustSEI[1])){ x$HC3RobustSEI[[cftname]]$pmap> -log10(0.05)} else NA })
    er = er[ !is.na(er)]
    blankimg[,,] = rowMeans(simplify2array(er), dims=3)
    allimgout$robustStatMapSEIpower[rowInd] = list(blankimg)
  } else {
    cftname = paste0('cft.p', value)
    # t type 1 error/power maps
    blankimg[,,] = 0 
    er = lapply(simdirs$results[simdirs$n==n], function(x){if(!is.na(x$tSEIp[1])){ x$tSEIp[[cftname]]$pmap> -log10(0.05)} else NA })
    er = er[ !is.na(er)]
    blankimg[,,] = rowMeans(simplify2array(er), dims=3)
    allimgout$tStatMapSEIpowerP[rowInd] = list(blankimg)
    # robust type 1 error/power maps
    er = lapply(simdirs$results[simdirs$n==n], function(x){if(!is.na(x$HC3RobustSEIp[1])){ x$HC3RobustSEIp[[cftname]]$pmap> -log10(0.05)} else NA })
    er = er[ !is.na(er)]
    blankimg[,,] = rowMeans(simplify2array(er), dims=3)
    allimgout$robustStatMapSEIpowerP[rowInd] = list(blankimg)
  }
  
}

#sapply(allimgout[ allimgout$value==0.25, 'robustStatmapMean'], function(x){x = x[,,30]; mask=mask[,,30]; quantile(x[ mask==1]) } )
#0.25^2*(simConfig$ns-ncol(design$X)) + design$df
save(allimgout, simConfig, template, file=file.path(dirname(simConfig$output), 'EST_reproducible_figs.rdata') )


## ----figure1, fig.width=11, fig.height=4.5------------------------------------
sval = 0.20
pval = 0.001
slices = 18:44

for(sval in simConfig$cfts.s){
  for(pval in simConfig$cfts.p){
    for(slice in slices){
      for(statmap in c('robustStatmapMean', 'tStatmapMean')){
        subimgout = allimgout[(allimgout$method=='p' & allimgout$value==pval) | (allimgout$method=='S' & allimgout$value==sval),]
        fname = file.path(dirname(simConfig$output), 'images', paste0(statmap, '_p', pval, '_s', sval, '_slice', slice, '.png') )
        dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
        png(filename = fname, height=4.5, width=11, units = 'in', res = 300)
        cex=1.5
        # graphical parameters
        fgcol = 'white'
        bgcol = 'black'
        par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
            cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
            mar=c(0.5,2.2,2,0), bty='l', oma=c(0.5,0,2,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
        layout(matrix(1:12, nrow=2, byrow=TRUE))
        for(rowInd in 1:nrow(subimgout)){
          n = subimgout[rowInd, 'n']
          method = subimgout[rowInd, 'method']
          value = subimgout[rowInd, 'value']
          
          design = pbj::getDesign(simConfig$form, simConfig$formred, data=simConfig$dat)
          if(method=='p'){
            chisqValue = qchisq(value, df=design$df, lower.tail=FALSE)
            othermethod = 'S'
            othervalue = round(sqrt(pmax((chisqValue - design$df)/n, 0)), 3)
          } else {
            #chisqValue = value^2*(n-ncol(design$X)) + design$df
            chisqValue = value^2*(n) + design$df
            othermethod = 'p'
            othervalue = round(pchisq(chisqValue, df=design$df, lower.tail = FALSE), 3)
          }
          if(othervalue==0){
            othervalue = '<0.001'
          } else {
            othervalue = paste0('=', othervalue)
          }
          blankimg = subimgout[[rowInd, statmap]]
          image(blankimg, template, thresh=chisqValue, index=slice, cex=cex*0.7)
          title(paste0('n=', n, ', ', othermethod, othervalue))
          if(rowInd == 1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
          if(rowInd == length(unique(subimgout$n))+1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
        }
        # if(statmap=='tStatmapMean'){
        #   mtext('T-statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # } else {
        #   mtext('Robust statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # }
        mtext('Target', side=3, outer = TRUE, cex=1*cex, font=2)
        dev.off()
      }
    }
  }
}



## ----figure2, fig.width=11, fig.height=4.5------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, round(ticks, 2), las=1, cex.axis=cex*0.7, font=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
# color thresholding for visualization and colorbar
threshs = c(0.05, 0.8)
cex=1.5

# graphical parameters
fgcol = 'white'
bgcol = 'black'

for(sval in simConfig$cfts.s){
  for(pval in simConfig$cfts.p){
    for(slice in slices){
      for(statmap in c("robustStatmapPower", "tStatmapPower")){
        subimgout = allimgout[(allimgout$method=='p' & allimgout$value==pval) | (allimgout$method=='S' & allimgout$value==sval),]
        fname = file.path(dirname(simConfig$output), 'images', paste0(statmap, '_p', pval, '_s', sval, '_slice', slice, '.png') )
        dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
        png(filename = fname, height=4.5, width=11, units = 'in', res = 300)
        layout(cbind(matrix(1:12, nrow=2, byrow=TRUE) %x% matrix(1, nrow=3, ncol=3) , 13))
        par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
            cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
            mar=c(0.5,2.2,2,0), bty='l', oma=c(0.5,0,2,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
        for(rowInd in 1:nrow(subimgout)){
          n = subimgout[rowInd, 'n']
          method = subimgout[rowInd, 'method']
          value = subimgout[rowInd, 'value']
          
          design = pbj::getDesign(simConfig$form, simConfig$formred, data=simConfig$dat)
          if(method=='p'){
            chisqValue = qchisq(value, df=design$df, lower.tail=FALSE)
            othermethod = 'S'
            othervalue = round(sqrt(pmax((chisqValue - design$df)/n, 0)), 3)
          } else {
            #chisqValue = value^2*(n-ncol(design$X)) + design$df
            chisqValue = value^2*(n) + design$df
            othermethod = 'p'
            othervalue = round(pchisq(chisqValue, df=design$df, lower.tail = FALSE), 3)
          }
          if(othervalue==0){
            othervalue = '<0.001'
          } else {
            othervalue = paste0('=', othervalue)
          }
          blankimg = subimgout[[rowInd, statmap]]
          image(blankimg, template, thresh=threshs, index=slice, cex=cex*0.7)
          title(paste0('n=', n, ', ', othermethod, othervalue))
          if(rowInd == 1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
          if(rowInd == length(unique(subimgout$n))+1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
        }
        # if(statmap=='tStatmapPower'){
        #   mtext('T-statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # } else {
        #   mtext('Robust statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # }
        mtext('Probability', side=3, outer = TRUE, cex=1*cex, font=2)
        par(mar=c(10,2.5,10,0.5))
        color.bar(pbj:::redyellow(64), min=threshs[1], max=threshs[2], nticks=4)
        dev.off()
      }
    }
  }
}


## ----figure3, fig.width=11, fig.height=4.5------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, round(ticks, 2), las=1, cex.axis=cex*0.7, font=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
# color thresholding for visualization and colorbar
threshs = c(0.05, 0.8)
cex=1.5

# graphical parameters
fgcol = 'white'
bgcol = 'black'

for(sval in simConfig$cfts.s){
  for(pval in simConfig$cfts.p){
    for(slice in slices){
      for(statmap in c("robustStatMapSEIpower", "tStatMapSEIpower", "robustStatMapSEIpowerP", "tStatMapSEIpowerP")){
        subimgout = allimgout[(allimgout$method=='p' & allimgout$value==pval) | (allimgout$method=='S' & allimgout$value==sval),]
        fname = file.path(dirname(simConfig$output), 'images', paste0(statmap, '_p', pval, '_s', sval, '_slice', slice, '.png') )
        dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
        png(filename = fname, height=4.5, width=11, units = 'in', res = 300)
        layout(cbind(matrix(1:12, nrow=2, byrow=TRUE) %x% matrix(1, nrow=3, ncol=3) , 13))
        par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
            cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
            mar=c(0.5,2.2,2,0), bty='l', oma=c(0.5,0,2,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
        for(rowInd in 1:nrow(subimgout)){
          n = subimgout[rowInd, 'n']
          method = subimgout[rowInd, 'method']
          value = subimgout[rowInd, 'value']
          
          design = pbj::getDesign(simConfig$form, simConfig$formred, data=simConfig$dat)
          if(method=='p'){
            chisqValue = qchisq(value, df=design$df, lower.tail=FALSE)
            othermethod = 'S'
            othervalue = round(sqrt(pmax((chisqValue - design$df)/n, 0)), 3)
          } else {
            #chisqValue = value^2*(n-ncol(design$X)) + design$df
            chisqValue = value^2*(n) + design$df
            othermethod = 'p'
            othervalue = round(pchisq(chisqValue, df=design$df, lower.tail = FALSE), 3)
          }
          if(othervalue==0){
            othervalue = '<0.001'
          } else {
            othervalue = paste0('=', othervalue)
          }
          blankimg = subimgout[[rowInd, statmap]]
          image(blankimg, template, thresh=threshs, index=slice, cex=cex*0.7)
          title(paste0('n=', n, ', ', othermethod, othervalue))
          if(rowInd == 1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
          if(rowInd == length(unique(subimgout$n))+1){ mtext(paste0(method, '=', value), side=2, cex=cex, font=2)}
        }
        # if(statmap=='tStatmapPower'){
        #   mtext('T-statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # } else {
        #   mtext('Robust statistic', side=3, outer = TRUE, cex=1*cex, font=2)
        # }
        mtext('Voxel-wise Type 1 error/Power', side=3, outer = TRUE, cex=1*cex, font=2)
        par(mar=c(10,2.5,10,0.5))
        color.bar(pbj:::redyellow(64), min=threshs[1], max=threshs[2], nticks=4)
        dev.off()
      }
    }
  }
}

