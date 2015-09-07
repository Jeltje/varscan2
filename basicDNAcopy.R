library(DNAcopy)

args <- commandArgs(TRUE)

# arguments are input and SD

# Alway use the same random seed for reproducible results
set.seed(0xcafe)	# cafe is a hex number

cn <- read.table(args[1],header=TRUE)
sd <- as.double(args[2])
CNA.object <-CNA(genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_stop, 
	data.type = 'logratio', sampleid = "sample")

smoothed.CNA.object <- smooth.CNA(CNA.object)

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=sd, verbose=1)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)

outfile <- paste(args[1], "SD", sd, "dnacopy.out", sep=".")

write.table(p.segment.smoothed.CNA.object[,1:6], file=outfile, quote=F, row.names=F, sep="\t")

detach(package:DNAcopy)
