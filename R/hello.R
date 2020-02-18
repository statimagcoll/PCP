

redisStart = function(){
  system('redis-server', wait=FALSE)
  library(rredis)
  redisConnect(nodelay = FALSE) #, password = "myredispw")
  # tests redis server
  #redisSet("test", runif(10))
  #redisGet("test")
}

# configures AWS using access key ID and secret access key provided in csv by user.
# Need to run on master. This sets the default access key and id.
# to get access key, click on name in top right of the management page, click security credentials, click access keys, click create access key.
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





# on fleet machines
# map drive where data are stored
# within docker run R CMD directing worker to the master node
# run docker image on AMI with job host argument pointing to master machine with -v to point to data available to AMI


# need to cancel spot fleet requests after running
aws ec2 cancel-spot-fleet-requests
