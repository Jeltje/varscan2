library(DNAcopy)

args <- commandArgs(TRUE)

# arguments are input, output basename, and chromosome

# read the full copycalled table, then select the chromosome we want
cn <- read.table(args[1],header=TRUE)
selectchrom <- args[3]
cn = cn[cn$chrom==selectchrom,]
CNA.object <-CNA(genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_stop, data.type = 'logratio', sampleid=c("Chromosome 9"))

smoothed.CNA.object <- smooth.CNA(CNA.object)

# setup output files
outfile <- paste(args[2],"dnacopy.out", sep=".")
outsd <- paste(args[2],"dnacopy.sd", sep=".")
outpng <- paste(args[2],"dnacopy.png", sep=".")
png(outpng, height=600, width=800)
par(mar=c(4,4,2,2))

# create iterable list for sdundo that goes from 4 to 0.5
sdlist<-seq(4, 0.5, -0.5)

for (i in sdlist){
#    write(paste("now at SD", i, sep=' '), stdout())
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=i, verbose=1)
    p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)
#    x <- c(i, length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]))
#    write(paste("At SD", i, x, "segments", sep=' '), file=outsd, append=TRUE)
#    write(x, file=outsd, append=TRUE)

# this tests how many output lines have $num.mark>=20. If there are fewer than 50, 
# it changes the sdundo and runs again
#    write(length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]), stdout())
    nrok <- length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]) 
    write(paste("At SD", i, nrok, "segments passed", sep=' '), file=outsd, append=TRUE)
    if(nrok >= 50){
	plot(segment.smoothed.CNA.object$data$maploc, segment.smoothed.CNA.object$data$Chromosome.9, pch=19, cex=0.25, cex.axis=1.25, cex.lab=1.5, col="cornflowerblue", ylim=c(-5,5), main="Chromosome 9", xlab="Position", ylab="Copy Number Change (log2)")
	segments(segment.smoothed.CNA.object$output$loc.start, segment.smoothed.CNA.object$output$seg.mean, segment.smoothed.CNA.object$output$loc.end, segment.smoothed.CNA.object$output$seg.mean, col="red", lwd=2)
	write.table(p.segment.smoothed.CNA.object, file=outfile, sep="\t")
	break
    }
}


detach(package:DNAcopy)
dev.off()
