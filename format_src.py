#!/usr/bin/env python
import numpy as np
import csv
import sys
import pylab

def read_list(filename):
	f=open(filename,'r')
	data = np.genfromtxt(f, skip_header=4)
	f.close()
	return data


# Read in output of imspect as I, Q, U files and output into a new complete file ready for RM code
print sys.argv[1]
if sys.argv[1] == '-h':
	print "Usage: format_src.py <srcname>  <rms>\n"
elif sys.argv[1]==None:
	print "Usage: format_src.py <srcname>  <rms>\n"
	exit
else:
	src = sys.argv[1]
	rms = float(sys.argv[2])
# File list names
stokespars = ('i','q','u')
il,ql, ul = ['%s.%s.log' % (src,a) for a in stokespars]
print "Reading %s %s %s " % (il, ql, ul)
idat = read_list(il)
qdat = read_list(ql)
udat = read_list(ul)
outfile="%s.dat" % src
f=open(outfile,'w')

for i in range(len(idat)):
	if idat[i,2]==0.0:
		print "Skipping channel %d" % i
	else:
		outstr = "%e  %e  %e  %e  %e  %e  %e\n" % (idat[i,1]*1.e9,idat[i,2],qdat[i,2],udat[i,2],rms,rms,rms)
#		print outstr
		f.write(outstr)

f.close()
