FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

# create a working directory and work from there
RUN mkdir /tmp/install
WORKDIR /tmp/install

RUN apt-get update && apt-get install -y \
	gcc \
	make \
	zlib1g-dev \
        git \
	wget \
	python-numpy \
	default-jre \
	r-base \
	bc

# DNAcopy version keeps changing so deprecated this:
# R and DNAcopy package (move to R library location)
#RUN apt-get install -y r-base
#RUN wget http://www.bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.40.0.tar.gz
# instead:
COPY ./DNAcopy_1.48.0.tar.gz ./ 
RUN R CMD INSTALL DNAcopy_1.48.0.tar.gz 

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
RUN tar -xvf samtools-0.1.18.tar.bz2
# the make command generates a lot of warnings, none of them relevant to the final samtools code, hence 2>/dev/null
RUN (cd samtools-0.1.18/ && make DFLAGS='-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=0' LIBCURSES='' 2>/dev/null && mv samtools /usr/local/bin)

# get varscan 
RUN wget  -O /usr/local/bin/VarScan.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar

# Move wrapper and helper scripts to same location
ADD ./run_varscan /usr/local/bin/
ADD ./separateArms.py /usr/local/bin/
ADD ./basicDNAcopy.R /usr/local/bin/
ADD ./meanLogRatioByChromosome.py /usr/local/bin/

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/install

ENTRYPOINT ["bash", "/usr/local/bin/run_varscan"]


