FROM ubuntu:14.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

RUN apt-get update && apt-get install -y \
        git \
	wget

# Create a new source directory
WORKDIR /home
RUN mkdir varscan-src
WORKDIR /home/varscan-src

# get varscan and run_varscan tools
RUN git clone https://github.com/dkoboldt/varscan.git
RUN git clone https://github.com/Jeltje/varscan2.git

# Move varscan.jar (allow for version updates) to /home directory
RUN rm varscan/VarScan.*source.jar
RUN mv varscan/VarScan.v2.*.jar /home/VarScan.jar

# move run_varscan scripts to /home directory
RUN mv varscan2/run_varscan /home
RUN mv varscan2/separateArms.py /home
RUN mv varscan2/basicDNAcopy.R /home
RUN mv varscan2/meanLogRatioByChromosome.py /home

# R and DNAcopy package (move to R library location)
RUN wget -O /home/DNAcopy_1.40.0.tar.gz http://www.bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.40.0.tar.gz

# Samtools 0.1.18 - note: 0.1.19 and 1.1 do NOT work, VarScan copynumber dies on the mpileup
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
RUN tar -xvf samtools-0.1.18.tar.bz2
# the make command generates a lot of warnings, none of them relevant to the final samtools code, hence 2>/dev/null
RUN (cd samtools-0.1.18/ && make DFLAGS='-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=0' LIBCURSES='' 2>/dev/null && mv samtools /home)

# Set working directory
WORKDIR /home
 
# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /home/varscan-src

