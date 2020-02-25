#!/bin/zsh
# This is run on the local machine (e.g. your laptop) to start an rstudio docker image
# start docker, in case you haven't already
#eval "$(docker-machine env default)"
# start rstudio with pbj docker instance. Port forwarding for rstudio and redis mapping dropbox drive
# 8787 is port for rstudio
# 6379 is port for redis
# runs in background
echo http://$(docker-machine ip default):8787
docker run -p 8787:8787 -p 6379:6379 -e ROOT=TRUE -e DISABLE_AUTH=true pbj:latest
# gets docker ip address
#ipadd=$(docker inspect --format '{{ .NetworkSettings.IPAddress }}' $(docker ps -a | head -2 | tail -1 | cut -d" " -f1) )

# ip address to connect to local rstudio instance
# From here: https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image


#!/bin/bash
# # Not run
# # connect to dropbox. Requires user interaction to copy web address into browser
# if [ ! -d ~/Dropbox ]; then
#   if [ ! -e ~/dropbox.tar.gz ]; then
#     wget -O dropbox.tar.gz "http://www.dropbox.com/download/?plat=lnx.x86_64"
#     tar -xvzf dropbox.tar.gz
#     rm -f dropbox.tar.gz
#   fi
#   .dropbox-dist/dropboxd
# fi
# chmod 600 ~/Dropbox/aws/pems/*pem
# # End not run
# sudo apt update
# sudo snap install docker
# sudo apt install -y libfuse2
# sudo apt install -y python3-pip
# pip3 install dbxfs
# export PATH=$PATH:/home/ubuntu/.local/bin
# echo "export PATH=$PATH:/home/ubuntu/.local/bin" >> ~/.bashrc
# mkdir ~/Dropbox
# dbxfs ~/Dropbox
# echo "dbxfs ~/Dropbox" >> ~/.bashrc
# cp /home/ubuntu/Dropbox/docker/dockerimages/pbj.tar ~/
# docker load -i ~/pbj.tar

# sudo apt update
# sudo snap install docker
# sudo apt install -y python3-pip
# wget https://github.com/dropbox/dbxcli/releases/download/v3.0.0/dbxcli-linux-amd64
# mkdir -p /home/ubuntu/.local/bin
# echo "export PATH=$PATH:/home/ubuntu/.local/bin" >> ~/.bashrc
# sudo chmod 711 dbxcli-linux-amd64
# mv dbxcli-linux-amd64 /home/ubuntu/.local/bin/dbxcli
# dbxcli get docker/dockerimages/pbj.tar ~/
# docker load -i pbj.tar
#
# mkdir ~/dropbox 2>/dev/null
# docker run -p 8787:8787 -p 6379:6379 -e ROOT=TRUE  -v "/home/ubuntu/dropbox:/home/rstudio/dropbox" -e DISABLE_AUTH=true pbj:latest
# docker run -p 8787:8787 -p 6379:6379 -e ROOT=TRUE  -e DISABLE_AUTH=true pbj:latest


# To increase memory available to docker open virtualbox and close machine. Then adjust memory and restart machine. Then open new terminal and start new docker image.
# This will be different on linux.

VBoxManage controlvm "default" poweroff
VBoxManage modifyvm "default" --memory 8192
VBoxManage startvm "default" --type headless


