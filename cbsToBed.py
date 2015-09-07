#!/usr/bin/python

import sys, os, re, getopt
usage = sys.argv[0]+""" <file>

Create bed file based on Varscan CBS output with
segment scores converted to colors (red for amplified, blue for deleted)

Option: 
	-n Don't create bed header
	-c <float> cutoff for amplification or deletion (default 0.25)

"""


# Main
# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "fdc:hbn")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)

makeBed = True
head = True
val = 0.25
for o, a  in opts:
#    if o == "-d":
#        doNotDelete = True
    if o == "-n":
        head = False
    if o == "-c":
        val = float(a)
    if o == "-h":
        print usage
        sys.exit()


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
    if fields[0] == 'ID':
	continue
    counter+=1
    id=('.').join(['mrg', str(counter)])
    qualifier = "neutral"
    rgb='0,0,0'
    if(float(fields[5]) < -val):
        rgb='0,0,255'	# blue
    elif(float(fields[5]) > val):
        rgb='255,0,0'	# red
    else:
        rgb='0,0,0'	# black
    chr = ("").join(["chr", fields[1]])
    outstring = ("\t").join([chr, fields[2], fields[3], id, '0', '.', fields[2], fields[2], rgb ])
    print outstring
f.close()


