#!/usr/bin/env python
import numpy as np
import os
#import astropy
import csv
from astropy import units as u
from astropy.coordinates import SkyCoord
import aplpy


def read_list(filename):
	with open(filename, 'r') as f:
		next(f) # skip header
		r = csv.reader(f, delimiter=" ",skipinitialspace=True)
		mydict = {rows[0]:rows[1:] for rows in r}
                #d = dict((row[0]), map(int,row[1:])) for row in r)
	return mydict

# Get the list of sources
f=open("src_pos.txt",'r')
sources = f.readlines()
f.close()

# File of RM components and sigma
fo=open('RC_pol.dat','wr')
fo.write("# SOURCE            RA              Dec      RM1      Err1P Err1M    RM2     Err2P Err2M      sigmaRM   \n")
basedir = 'data/'
rm1,rm1errp,rm1errm = (0.0,0.0,0.0)
for s in sources:
	source = s.split()[0]
	ra = s.split()[1]
	dec = s.split()[2]
# Initialise some values
	rm1,rm1errp,rm1errm = (0.0,0.0,0.0)
	rm2,rm2errp,rm2errm = (0.0,0.0,0.0)
	sigRM,sigRMerrp,sigRMerrm = (0.0,0.0,0.0)
# Open the result file
	result_file = basedir+source+"_QUfit.result"
	if os.path.exists(result_file):
		f=open(result_file)
		print "Opening %s" %result_file
	else: 	
		f=open(basedir+"1comp/"+source+"_QUfit.result")
		print "Opening "+basedir+"1comp/"+source+"_QUfit.result"
	for line in f.readlines():	
		if line.startswith("RM1_"):
			rm1 = float(line.split()[2])
			rm1errp = float(line.split()[3][1:-1])
			rm1errm = float(line.split()[4][1:-1])
		elif line.startswith("RM_"):
			rm1 = float(line.split()[2])
			rm1errp = float(line.split()[3][1:-1])
			rm1errm = float(line.split()[4][1:-1])
		elif line.startswith("RM2_"):
			rm2 = float(line.split()[2])
			rm2errp = float(line.split()[3][1:-1])
			rm2errm = float(line.split()[4][1:-1])
		elif line.startswith("sigmaRM_"):
			sigRM = float(line.split()[2])
			sigRMerrp = float(line.split()[3][1:-1])
			sigRMerrm = float(line.split()[4][1:-1])

	outstr = "%s  %s  %s   %.2f     %.2f %.2f     %.2f     %.2f %.2f       %.2f     %.2f %.2f\n" %(source,ra,dec,rm1,rm1errp,rm1errm,rm2,rm2errp,rm2errm,sigRM,sigRMerrp,sigRMerrm)
	fo.write(outstr)
	f.close()


# Now open that file and read in the data
mydict=read_list("RC_pol.dat")
rastr = []
decstr = []
rm1 = []
rm1err = []
rm2 = []
sigmaRM = []
exclude=[]
lat=[]
lon=[]
markersize=[]
sources = []
# Plot the RMs on the data
input_cube='/Users/naomi/Data/GC/GC.hi.sub.fits'
input_cube='/Users/naomi/Data/GASS/Galcen/GASS_GC.sub.fits'
fig=aplpy.FITSFigure(input_cube,slices=[33],dimensions=[0,1])
fig.show_grayscale(vmin=15,vmax=110)
fig.show_colorbar()
fig.add_label(0.5, 1.03, 'Parkes data',layer='title',relative=True)
fig.recenter(3,1,width=17.,height=15.)


# Work out the sources to exclude based on the RM error (currently set to greater than 10)
for src in mydict.keys():
	if float(mydict[src][3])>=10.0:
		exclude.append(src)

for src in mydict.keys():
	if not src in exclude:		
		sources.append(src)
		rastr.append(mydict[src][0])
		decstr.append(mydict[src][1])
		rm1.append(float(mydict[src][2]))
		rm2.append(float(mydict[src][5]))
		sigmaRM.append(float(mydict[src][8]))
		rm1err.append(float(mydict[src][3]))
		coord=SkyCoord(ra=rastr, dec=decstr,unit=(u.hourangle,u.degree),frame='icrs')
lon=coord.galactic.l.degree
lat=coord.galactic.b.degree

for i in range(len(lon)):
	print sources[i],lon[i],lat[i],rm1[i],rm2[i]
	if rm1[i]>0.0:
		print "Red marker "
		fig.show_markers(lon[i],lat[i], edgecolor='red', facecolor='red',marker='+', s=10.*abs(rm1[i]))
	elif rm1[i]<0.0:
		print "Blue marker "
		fig.show_markers(lon[i],lat[i], edgecolor='blue', facecolor='blue',marker='+', s=10.*abs(rm1[i]))

# Now plot the second RM component
	if rm2[i]>0.0:
		fig.show_markers(lon[i],lat[i], edgecolor='red',marker='o', s=10.*abs(rm2[i]))
	elif rm2[i]<0.0:
		fig.show_markers(lon[i],lat[i], edgecolor='blue', marker='o', s=10.*abs(rm2[i]))


