#!/inside/home/cline/bin/python

import sys, os, re, argparse

parser = argparse.ArgumentParser(description="Divide input copycaller file into p and q segments using an input bed file with centromere coordinates")
parser.add_argument('input', type=str,help="Tab separated copycaller file")
parser.add_argument('centro', type=str,help="Centromere bed file")
# optional argument
#parser.add_argument('-t', dest='targetfile', type=str, default=False,
#        help="Input file with target genes and scores")
# optional flag
#parser.add_argument('-d', '--debug', help="Optional debugging output", action='store_true')

def myfunction(line):
    """DESCRIBE ME"""
    line = line[:-1]

class cent(object):
    """centromere objects by chromosome"""
    def __init__(self, line):
        [id, start, end] = line.split("\t")
	self.name = id
        self.p = int(start)
        self.q = int(end)

def getChrom(id, centList, cur):
    if cur and id == cur.name:
	return cur
    for i in centList:
	if id == i.name:
	    return i
    return False

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Main

centros=[]
f = open(args.centro,'r')
for line in f:
    line = line.strip()
    centobj = cent(line)
    centros.append(centobj)
f.close

curchrom = False
f = open(args.input,'r')
for line in f:
    line = line.strip()
    fields = line.split("\t")
    if fields[0] == 'chrom':
	print line
	continue
    id = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    curchrom = getChrom(id, centros, curchrom)
#    print >>sys.stderr, id, curchrom
#    if (not curchrom) or (not curchrom.name == fields[0]):
#	curchrom = False
#	for c in centros:
#	    if c.name == id:
#		curchrom = c
#		print >>sys.stderr, "setting", curchrom
#		break
    # centromere free chromosomes
    if not curchrom:
	print line
    # print only if segment does not overlap centromere
    elif end < curchrom.p:
	fields[0] = ('.').join([fields[0], 'p'])
	print ("\t").join(fields)
    elif start > curchrom.q:
	fields[0] = ('.').join([fields[0], 'q'])
	print ("\t").join(fields)
f.close
