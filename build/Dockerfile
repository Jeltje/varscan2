FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
	gcc \
	make \
	zlib1g-dev \
        git \
	wget

# Create a new source directory
WORKDIR /home
RUN mkdir exec

# get varscan 
RUN wget  -O /home/exec/VarScan.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar

# R and DNAcopy package (move to R library location)
# NOTE: this is NOT a stable way of doing this, bioconductor updates this package regularly
RUN wget -O /home/exec/DNAcopy_1.42.0.tar.gz http://www.bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.46.0.tar.gz

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
RUN tar -xvf samtools-0.1.18.tar.bz2
# the make command generates a lot of warnings, none of them relevant to the final samtools code, hence 2>/dev/null
WORKDIR /home/samtools-0.1.18
RUN make DFLAGS='-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=0' LIBCURSES='' 2>/dev/null 
RUN mv samtools /home/exec

# change ownership of outputfiles so they can be deleted
RUN chmod -R 777 /home/exec

# Set working directory
WORKDIR /home
 
# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /home/samtools*

