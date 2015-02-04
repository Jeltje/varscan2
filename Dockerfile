FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

# create a working directory and work from there
RUN mkdir /tmp/install
WORKDIR /tmp/install

# This runs apt-get without confirmations
# but so does -y (assume yes), so this might not be necessary
# ENV DEBIAN_FRONTEND=noninteractive

# make sure the package repository is up to date
RUN apt-get -qq update
RUN apt-get install -y wget
RUN apt-get install -y bc

# R and DNAcopy package (move to R library location)
RUN apt-get install -y r-base
RUN wget http://www.bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.40.0.tar.gz
RUN R CMD INSTALL DNAcopy_1.40.0.tar.gz 

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
RUN tar -xvf samtools-0.1.18.tar.bz2
RUN (cd samtools-1.1/ && make && make prefix='/usr/' install)

#java
RUN apt-get install -y default-jre
RUN wget -O /opt/VarScan.v2.3.7.jar http://sourceforge.net/projects/varscan/files/VarScan.v2.3.7.jar/download

# Varscan scripts
# for now let's mount run_varscan.sh from outside the image so we can tweak it
RUN wget -O /usr/local/bin/iterDNAcopy.R http://hgwdev.cse.ucsc.edu/~jeltje/iterDNAcopy.R
RUN chmod a+x /usr/local/bin/iterDNAcopy.R

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* 

# TODO Set default container command
# ENTRYPOINT usr/local/bin/run_varscan
