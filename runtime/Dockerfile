FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
        default-jre \
	r-base \
        bc


# R and DNAcopy package (move to R library location)
ADD DNAcopy_1.40.0.tar.gz
RUN R CMD INSTALL DNAcopy_1.40.0.tar.gz 

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
ADD samtools /usr/local/bin

# Copy over jar and place in /opt/cgl-docker-lib
RUN mkdir /opt/cgl-docker-lib
ADD VarScan.jar  /opt/cgl-docker-lib/

# Move wrapper and helper scripts to same location
ADD run_varscan /opt/cgl-docker-lib
ADD separateArms.py /opt/cgl-docker-lib
ADD basicDNAcopy.R /opt/cgl-docker-lib
ADD meanLogRatioByChromosome.py /opt/cgl-docker-lib

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/cgl-docker-lib/run_varscan"]

