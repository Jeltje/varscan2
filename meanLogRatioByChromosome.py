#!/usr/bin/python

import sys, os, re, argparse
from numpy import *

parser = argparse.ArgumentParser(description="Calculates log average per chromosome")
parser.add_argument('cpcalled', type=str,help="copyCalled output")

class Chrom(object):
    """Holds fragment scores by chromosomes """
    def __init__(self, chrom, score, endpos):
        self.frags = []
	self.chrom = chrom
	self.add(chrom, score, endpos)
    def add(self, chrom, score, endpos):
	if self.chrom == chrom:
	    self.frags.append(score)
	    self.lastfrag = endpos
	    return True
	else:
	    return False
    def stats(self):
	self.mean = mean(self.frags)
	self.std = std(self.frags)


if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Main
chromTable = []	# holds chrom objects
curChrom = Chrom('empty', 0, 0)
f = open(args.cpcalled,'r')
for line in f:
    line = line.strip()
    fields = line.split("\t")
    chr = fields[0]
    if chr == "chrom":
	continue
    score = float(fields[6])
    endpos = int(fields[2])
    if not (curChrom.add(chr, score, endpos)):
	curChrom = Chrom(chr, score, endpos)
	chromTable.append(curChrom)
f.close()

if len(chromTable) < 3:
    print >>sys.stderr, "ERROR: Please enter whole genome file"
    sys.exit(1)

means = []
for chr in chromTable:
    if chr.lastfrag < 60000000:	# skip MT, GL, Y
#    if len(chr.frags) < 3000:   # skip MT, GL, Y
	continue
    chr.stats()
    means.append(chr.mean)

med = len(means)/2   # this rounds, which is perfect
means.sort()
print "%.2f" % mean([means[med-1], means[med], means[med+1]])

