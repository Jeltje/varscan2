#!/usr/bin/python

import sys, os, re, getopt
usage = sys.argv[0]+""" <file>

Create bed file based on mergeSegments output with
amplified and deleted codes converted to colors

Option: 
	-n Don't create bed header

"""


# Main
# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "fdchbn")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)

makeBed = True
head = True
for o, a  in opts:
#    if o == "-d":
#        doNotDelete = True
    if o == "-n":
        head = False
    if o == "-h":
        print usage
        sys.exit()

#some = 0.249
#very = 1

if len(args) != 1:
    sys.exit(usage)

# Run program

f = open(args[0],'r')
counter = 0

if makeBed and head:
    print 'track name=%s description="%s" itemRgb="On"'% (args[0], args[0])
for line in f:
    line = line.strip()
    fields = line.split("\t")
    if fields[0] == 'chrom':
	continue
    counter+=1
    id=('.').join(['mrg', str(counter)])
    qualifier = "neutral"
    rgb='0,0,0'
    if(fields[7] == 'deletion'):
        rgb='0,0,255'	# blue
    elif(fields[7] == 'amplification'):
        rgb='255,0,0'	# red
    elif(fields[7] == 'neutral'):
        rgb='0,0,0'	# black
    else:
	print >>sys.stderr, "don't understand", fields[7]

    if makeBed:
        outstring = ("\t").join([fields[0], fields[1], fields[2], id, '0', '.', fields[1], fields[1], rgb ])
    else:
        outstring = ("\t").join([fields[0], fields[1], fields[2], fields[4], qualifier])
    print outstring
f.close()


