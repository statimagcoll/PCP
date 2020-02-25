
## ----configureAWS--------------------------------------------------------
 dbxcli = function(cmd, ..., ncores=1){
    result = parallel::mclapply(paste('dbxcli', cmd, do.call(paste, list(...) ) ), system, mc.cores=ncores )
  }

configureAWS = function(keycsv, profile='[default]', region='us-east-2', output='json'){
  dir.create('~/.aws', showWarnings = FALSE)
  if(!is.null(keycsv)){
    suppressWarnings(key <- read.table(keycsv, stringsAsFactors = FALSE, sep=',', header = FALSE))
    id = gsub('.*=', '', key[1,1])
    key = gsub('.*=', '', key[2,1])

    fileConn<-file("~/.aws/credentials")
    writeLines(paste0(c('', 'aws_access_key_id=', 'aws_secret_access_key='), c(profile, id , key) ), fileConn)
    close(fileConn)

    fileConn<-file("~/.aws/config")
    writeLines(paste0(c('', 'region=', 'output='), c(profile, region, output)), fileConn)
    close(fileConn)
  }
}

dbxcli('get',  'aws/aws_keys/vandeks.csv',  '~/')
awskey = '~/vandeks.csv'
configureAWS(awskey)


## ---- redisVariables-----------------------------------------------------
compute.config = list(
  queuename = 'pbj',
  host = gsub('^.*://|:.*$', '', Sys.getenv('RSTUDIO_HTTP_REFERER')),
  # password for this redis session
  password = paste(sample(c(0:9,letters,LETTERS), 16),collapse=""),
  ncores = 3
)
set.seed(666)


## ---- startRedisServer---------------------------------------------------
# default redis config template
# only option to change is password
# doRedis does not appear to support non protected mode
dbxcli('get', 'redis/redisTemplate.conf', '~/')
redisConfTemplate = '~/redisTemplate.conf'
redisConf = '~/redis.conf'
# insert password into file
system(paste0('sed s/###PASSWORD###/', compute.config$password, '/ ', redisConfTemplate, ' > ', redisConf))
system(paste0('redis-server ', redisConf), wait=FALSE)
library(doRedis)
registerDoRedis(compute.config$queuename, nodelay=FALSE, password=compute.config$password)


## ---- dockerHostTest, eval=FALSE-----------------------------------------
## library(doRedis)
## startLocalWorkers(n=compute.config$ncores, queue=compute.config$queuename, host=compute.config$host, password = compute.config$password,
##                   nodelay=FALSE, linger=60*10, timeout=60*10)


## ---- dockerTest, eval=FALSE, echo=FALSE---------------------------------
## test = foreach(i=1:10, .combine=c) %dopar% {
##   control = rnorm(200)
##   exp = rnorm(200)
##   t.test(control, exp)$p.value
## }


## ------------------------------------------------------------------------
### LIBRARIES ###
library(RNifti)
library(parallel)
library(splines)
library(mmand)
library(fslr)
library(progress)
library(abind)
library(pbj)


### LOAD IN DATA FROM DROPBOX ###
dbimagedir = 'pbj/data/abide/neuroimaging/cpac/alff'
dbdatafile = 'pbj/data/abide/demographic/n1035_phenotypic_20190509.rds'


# load in data and get directories
dbxcli('get', dbdatafile, '~/data.rds')
dat = readRDS('~/data.rds')
# fake covariates for later use if needed
dat$fake_covariate = rnorm(nrow(dat))
dat$fake_group = factor(sample(1:4, nrow(dat), replace=TRUE))

# same with imaging data
datadir = '~/data/images'
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dat$imgname = paste(dat$file_id, 'alff.nii.gz', sep='_')
dat$images = file.path(datadir, dat$imgname)
dbimages = file.path(dbimagedir, dat$imgname)
# download the images from dropbox
dbxcli('get', dbimages, dat$images, ncores=compute.config$ncores)
maskfile = '~/mask.nii.gz'
dbmaskfile = 'pbj/data/abide/neuroimaging/MNI152_T1_3mm.nii.gz'
dbxcli('get', dbmaskfile, maskfile)


### SIMULATION PARAMETERS ###
sim.config = list(
  robust = TRUE, # use robust variance estimator?
  # vector of sample sizes to simulate
  ns = c(50),
  # number of simulations to run, or vector of simulation   numbers to rerun
  nsim=100,
  # number of bootstraps
  nboot = 500,
  # number of permutations
  nperm = 500,
  # cluster forming thresholds
  cfts.s = c(0.1, 0.25, 0.4),
  cfts.p = c(0.05, 0.01, 0.005, 0.001),
  
  # radius for spheres of signal.
  rs=c(8),
  
  #### MODEL FORMULAS FOR SIMULATIONS ####
  formres = as.formula( paste0(" ~ dx_group + sex + func_mean_fd + ns(age_at_scan, df=10)" )),
  # need age_at_scan in both models for testing nonlinear functions
  form = as.formula(paste0(" ~ sex + func_mean_fd + age_at_scan + ns(fake_covariate, df=4)" )),
  formred = as.formula(paste0(" ~ sex + func_mean_fd + age_at_scan + fake_covariate")),
  # variance estimator variable for deweighting
  varvar = c("func_mean_fd")
)
# use betas = 0 for global null
# parameters = betas * sd(y)/sd(x).
sim.config$betas = rep(0, length(sim.config$rs))



### SETUP THE SIMULATION ANALYSIS ###
# subsets dataset to all people who have the variables
dat = dat[apply(!is.na(dat[ ,c(all.vars(as.formula(sim.config$formres)), sim.config$varvar)]), 1, all), ]
dat$rfiles = file.path('~/data/res_images', basename(dat$images))


pbj::residualizeImages(files=dat$images, dat=dat, mask=maskfile, form=sim.config$formres, outfiles=dat$rfiles, mc.cores=compute.config$ncores)

# output directories. Probably won't need these with doRedis
#name=paste(dataset, paste(gsub(" +", "_", gsub(" +$|^ +", "", gsub("[[:punct:]]", "", form))),
#                          gsub(" +", "_", gsub(" +$|^ +", "", gsub("[[:punct:]]", "", formred))), sep='-'), sep='_')
#form = as.formula(form)
#formred = as.formula(formred)











## ---- setupAWSjson-------------------------------------------------------
# output name of the json config file -- THIS FILE WILL BE MADE OR OVERWRITTEN
createLaunchJSON = '/home/rstudio/dropbox/aws/ec2/create_launch_template_jsons/pbjWorkers.json'
# get host ip address
createLaunchTemplateJSON = '/home/rstudio/dropbox/aws/ec2/create_launch_template_jsons/pbjWorkersTemplate.json'
#readLines(templatejson)

# get commands for UserData field
userdata = tempfile()
userdatacmds = paste("#!/bin/bash; R --file=~/home/ubuntu/start_redis_workers.R --args", host, queuename, collapse=' ' )
fileConn = file(tmpfile)
writeLines(userdata, fileConn)
close(fileConn)
#userdata = system(paste('base64', tmpfile), intern=TRUE)
# the host argument location is identified by the string ###HOST###
system(paste0('sed s+###USERDATA###+',tmpfile, '+ ', createLaunchTemplateJSON, ' > ', createLaunchJSON))


## ---- createLaunchTemplate-----------------------------------------------
# create spot instances
aws ec2 run-instances --image-id ami-09b4d361c6a48c1d5 --count 2 --instance-type t2.micro --key-name redisServer --user-data file:///tmp/Rtmpm5qHyi/file10a178470db --instance-market-options file:///home/rstudio/dropbox/aws/ec2/market_type/specification.json
# copy pem file over to spot instances

# update and install sshfs
# map dropbox drive from spot instances


## ---- createFleet, eval=FALSE--------------------------------------------
## # createFleet json config file
## #createFleetJSON = '/home/rstudio/dropbox/aws/ec2/create_fleet_jsons/pbjWorker.json'
## #fleetJSON = system(paste0('aws ec2 create-fleet --cli-input-json file://', createFleetJSON), intern=TRUE )
## #requestSpotInstancesJSON = '/home/rstudio/dropbox/aws/ec2/request_spot_instances/request_little_spot_instances.json'
## #spotsJSON = system(paste0('aws ec2 request-spot-instances --cli-input-json file://', requestSpotInstancesJSON), intern=TRUE )
## #aws ec2 request-spot-instances --instance-count 2 --type "one-time" --launch-specification file:///home/rstudio/dropbox/aws/ec2/request_spot_instances/specification.json
## # run for loop
## # close fleet

