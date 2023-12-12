# TB Hackday RStudio Server Launch Instructions

## Overall steps to launch
1. Launch AWS EC2 instances using the AWS console
2. Login to each instance via SSH
3. Pull and run docker container in the background
4. Attach to the docker container and clone the hackday github repos
5. Create the user accounts for the hackday (all accounts on all instances to flexibly balance load)
6. Download data from figshare
7. Test user access to RStudio. Post IP addresses to slack.
8. Detach from container and logout from EC2 instance


## Launch AWS EC2 instances
 - Login to the AWS console using 12 digit IAM key, username and password: AWS Console (https://aws.amazon.com/console/)
 - Check that you are in the `us-west-2` region
 - Launch existing instance that is pre-configured. Launch multiple to accomodate more users.
 - Save list of instance IP addresses to post
 - See below for details of how to configure a new instance

## Login to EC2 instance
 - Ensure that the `ssh-agent` has been started in current session (`eval $(ssh-agent -s)`)
 - Ensure that the correct PEM file has been added (e.g., `ssh-add ~/.ssh/agartlan-west-2.pem`)
 - Login with default username: `ssh ec2-user@[INSTANCE-IP]`

## Launch RStudio `docker` image
Pull the latest container from dockerhub.com and launch. Port is specified as `-p [host]:[container]`. Note that running docker container in interactive mode seems to eventually lead to an unstable RStudio session. Use `-d` instead of `-it` and detach once container is configured.

```
docker pull afioregartland/hackday-rstudio
docker run -d -p 8787:8787 -e PASSWORD=dude97 --name hackday_rstudio afioregartland/hackday-rstudio:latest
docker exec -it hackday_rstudio /bin/bash
```

## Setup the RStudio container

Login to the container and down load the latest code repos, create system users for RStudio and download the data. Order is important: create users then copy PAM profile. Restarting the RStudio service may help if users are not recognized. See PAM in this doc for help: [RStudio Server Admin Guide](https://s3.amazonaws.com/rstudio-server/rstudio-server-pro-0.98.1074-admin-guide.pdf). Users may need to chwd("/home") to see all the files.

```
cd /home
git clone https://github.com/FredHutch/seatrac-hackday-2023.git
git clone https://github.com/agartland/hackday-rstudio.git

python3 hackday-rstudio/create_users.py seatrac-hackday-2023/roster.users.csv
cat /etc/passwd
cp /etc/pam.d/login /etc/pam.d/rstudio

mkdir /home/bigdata
wget https://figshare.com/ndownloader/articles/24425053/versions/3 -O /home/bigdata/bigdata.zip
cd /home/bigdata
unzip bigdata.zip
```

It may be helpful to stop and restart the RStudio server at some point.
```
service --status-all
service rstudio-server restart
```

## Additional configuration notes for the AWS EC2 instance
Use a launch template for quickly launching a new instance. The template `RStudio-Server` currently has the following details.
 - Instance types: `c5.12xlarge` has 48 vCPUs, 96 GiB of RAM and costs $2.04/hr (Dec-2023)
 - AMI: AWS Linux 2023
 - Security group: `launch-wizard-1` currently has several inbound ports open for RStudio or jupyter notebooks including 8989 for RStudio default
 - SSH keys: Use existing key-pair and make sure PEM file is added to the ssh-user before trying to SSH
 - Storage should be sufficient for data by default, given how much memory we are procuring, but volume could be expanded
 - User data: This script is run when instance starts
 
 ```bash
 #! /bin/sh
sudo yum update -y
sudo yum -y install docker
sudo service docker start
sudo usermod -a -G docker ec2-user
chkconfig docker on
docker pull afioregartland/hackday-rstudio
```

## Other useful `docker` bits
 - `docker build --tag afioregartland/hackday-rstudio:latest https://github.com/agartland/hackday-rstudio.git`
 - `docker images`
 - `docker tag [IMAGE_ID] afioregartland/hackday-rstudio:latest`
 - `docker run -d -p 8787:8787 -e PASSWORD=dude97 --name hackday_rstudio afioregartland/hackday-rstudio:latest`
 - `docker exec -it hackday_rstudio /bin/bash`
 - `docker push afioregartland/hackday-rstudio:latest`

## Unit test for `lme4` installation

Both `Seurat` and `lme4` depend on `Matrix` and can be built against the latest version `1.6-4`. However, this was not the default version that came with the `rocker/tidyverse` container and so the Dockerfile was modified to remove the older version (`1.5-3`) and install the latest versions of these tools from source. This test ensures that `lme4` is working.

```
library(lme4)
library(palmerpenguins)

dat <- penguins
lme4::lmer(bill_length_mm ~ year + (1|species), data=dat)
```
