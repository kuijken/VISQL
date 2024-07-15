import numpy as np
import astropy.io.fits as pf
from sys import argv
from datetime import datetime
from astropy.coordinates import get_sun
from astropy.time import Time

# use time to get solar position vector s
# use target position vector t to work out SAA:
#     cos(SAA) = s.t
# then define vector a that points away from sun in focal plane:
# a = (s x t) x t    (normalised to unit length)

def radec2vec(radec):  # turn ra,dec (deg) into a unit vector, X:ra=0, Z:dec=90
    ra,dec=radec[0]*np.pi/180,radec[1]*np.pi/180
    xyz=np.zeros((3))
    xyz[0]=np.cos(ra)
    xyz[1]=np.sin(ra)
    xyz  *=np.cos(dec)
    xyz[2]=np.sin(dec)
    return xyz

def normvec(xyz):
    return xyz/(xyz**2).sum()**0.5

def getEuclidAngles_old(sunradec,targradec,pa):
    '''
    sun is (RA,DEC) of the sun
    targ is (RA,DEC) of the centre of VIS
    pa is the PA of VIS as returned by Astrometry.net
    fpy is vector in tel FP that points away from sun
    '''
    sunxyz=radec2vec(sunradec)
    targxyz=radec2vec(targradec)
    cosSAA=np.dot(sunxyz,targxyz)
    SAA=np.arccos(cosSAA)*180/np.pi
    fpy=normvec(np.cross(np.cross(sunxyz,targxyz),targxyz))
    fpE=normvec(np.cross(targxyz,[0,0,1]))  # points East on image
    alpha0=np.arccos(np.dot(fpy,fpE))*180/np.pi
    return SAA,alpha0

def getEuclidAngles(sunradec,targradec,pa):
    '''
    sun is (RA,DEC) of the sun
    targ is (RA,DEC) of the centre of VIS
    pa is the PA of VIS as returned by Astrometry.net
    fpy is vector in tel FP that points away from sun
    '''
    sunxyz=radec2vec(sunradec)
    targxyz=radec2vec(targradec)
    cosSAA=np.dot(sunxyz,targxyz)
    SAA=np.arccos(cosSAA)*180/np.pi
    fpy=normvec(np.cross(np.cross(sunxyz,targxyz),targxyz))
    fpE=normvec(np.cross(targxyz,[0,0,1]))  # points East on image
    fpN=normvec([0,0,1]-np.dot(targxyz,[0,0,1])*targxyz)
    #alpha0=np.arccos(np.dot(fpy,fpE))*180/np.pi
    alpha0=np.arctan2(np.dot(fpy,fpN),np.dot(fpy,fpE))*180/np.pi
    return SAA,alpha0

obsname = argv[1]     # png or fits file name
if obsname[-4]=='.': obsname=obsname[:-4]
if obsname[-5]=='.': obsname=obsname[:-5]

# extract obstime
yy=obsname[2:6]
mm=obsname[6:8]
dd=obsname[8:10]
h=obsname[11:13]
m=obsname[13:15]
s=obsname[15:17]
datestring=yy+"-"+mm+"-"+dd+"T"+h+":"+m+":"+s
dateobs=Time(datestring)
sun=get_sun(dateobs)
sunradec=np.array([sun.ra.degree,sun.dec.degree])

# extract target and PA
wcsname='../SCIENCE/'+obsname+'.wcs'
try:
    wfile=pf.open(wcsname)
    wh=wfile[0].header
    ratarg=wh['CRVAL1']
    detarg=wh['CRVAL2']
    cd1_1= wh['CD1_1']
    cd1_2= wh['CD1_2']
    cd2_1= wh['CD2_1']
    cd2_2= wh['CD2_2']
    targradec=[ratarg,detarg]
    pa=(np.arctan2(cd1_2,cd1_1)-np.arctan2(cd2_1,cd2_2))/2 * 180/np.pi
#    print('Target:',ratarg,detarg,pa)
    wfile.close()
except:
    print(wcsname, 'not found.')
    exit()

saa,alpha0=getEuclidAngles(sunradec,targradec,pa)
print(saa,(180-alpha0+pa) % 360 - 180,pa,ratarg,detarg)

#saa,alpha0=getEuclidAngles_old(sunradec,targradec,pa)
#print(saa,alpha0+pa,pa,ratarg,detarg)


'''
Looks like the alpha=alpha0+pa where
alpha = spacecraft roll angle
alpha0= angle from N to line in s/c FP pointing away from sun
pa    = orientation of image on sky from astrometry.net
'''
