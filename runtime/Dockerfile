FROM ubuntu:14.04

# Note: FROM java and FROM r-base work too but take much longer apt-get updating

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	python-numpy \
        default-jre \
	r-base \
        bc


# R and DNAcopy package (move to R library location)
# Note: cannot use ADD because this untars the file
COPY exec/DNAcopy_1.42.0.tar.gz /tmp/
RUN R CMD INSTALL /tmp/DNAcopy_1.42.0.tar.gz 

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
ADD exec/samtools /usr/local/bin/

# Copy over jar and place in /opt/cgl-docker-lib
RUN mkdir /opt/cgl-docker-lib
ADD exec/VarScan.jar  /opt/cgl-docker-lib/

# Move wrapper and helper scripts to same location
ADD ./run_varscan /opt/cgl-docker-lib/
ADD ./separateArms.py /opt/cgl-docker-lib/
ADD ./basicDNAcopy.R /opt/cgl-docker-lib/
ADD ./meanLogRatioByChromosome.py /opt/cgl-docker-lib/

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/DNAcopy_1.42.0.tar.gz

ENTRYPOINT ["bash", "/opt/cgl-docker-lib/run_varscan"]

