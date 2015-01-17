library(DNAcopy)

args <- commandArgs(TRUE)

#regions <- read.table("input/dnacopy.input")
#CNA.object <- CNA(regions$V6, regions$V1, regions$V2, data.type="logratio", sampleid=c("Chromosome 9"))

write(args[1], stdout())
png("output/dnacopy.jpg", height=600, width=800)
par(mar=c(4,4,2,2))

#Was
# from help(CNA):
#       CNA(genomdat, chrom, maploc, data.type=c("logratio","binary"),
#                      sampleid=NULL, presorted = FALSE)
# therefore

regions <- read.table(args[1])	
CNA.object <-CNA(genomdat = regions[,6], chrom = regions[,1], maploc = regions[,2], data.type = 'logratio')

# somehow these have to both be changed - you can't do $V6 when reading args[1]


write.table(CNA.object, file="output/CNA.out", sep="\t")
write(args[1], stdout())


# note that sampleid and chrom are both getting assigned and are likely the reason for the odd output fields

smoothed.CNA.object <- smooth.CNA(CNA.object)

# create iterable list for sdundo that goes from 4 to 0.5
sdlist<-seq(4, 0.5, -0.5)

for (i in sdlist){
    write(i, stdout())
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=i, verbose=1)
    p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)
    x <- c(i, length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]))
    write(x, file="output/dnacopy.segments.sd", append=TRUE)

# this tests how many output lines have $num.mark>=20. If there are fewer than 50, 
# it changes the sdundo and runs again
    write(length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]), stdout())
    if(length(p.segment.smoothed.CNA.object$pval[p.segment.smoothed.CNA.object$num.mark>=20]) >= 50){
	plot(segment.smoothed.CNA.object$data$maploc, segment.smoothed.CNA.object$data$Chromosome.9, pch=19, cex=0.25, cex.axis=1.25, cex.lab=1.5, col="cornflowerblue", ylim=c(-5,5), main="Chromosome 9", xlab="Position", ylab="Copy Number Change (log2)")
	segments(segment.smoothed.CNA.object$output$loc.start, segment.smoothed.CNA.object$output$seg.mean, segment.smoothed.CNA.object$output$loc.end, segment.smoothed.CNA.object$output$seg.mean, col="red", lwd=2)
#	write.table(p.segment.smoothed.CNA.object, file="output/test.out")
	write.table(p.segment.smoothed.CNA.object, file="output/test.out", sep="\t")
	break
    }
    write("continuing", stdout())
}


detach(package:DNAcopy)
dev.off()
