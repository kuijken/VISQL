# read fits header and write a header that should be very close to what astrometry.net makes for a thumbnail.

#import astropy.wcs as wcs
import astropy.io.fits as pf
from sys import argv
#from glob import glob
from os import path
from numpy import pi,cos,sin

fitsfiles=argv[1:]
wcsfiles=[x[:-4]+'wcs' for x in fitsfiles]
wcsfilestbd=[x for x in wcsfiles if not path.exists(x)]
deg=pi/180
for tbd in wcsfilestbd:
  try:
    oldh=pf.open(tbd[:-3]+'fits')[0].header
    ra  = oldh['RA']
    dec = oldh['DEC']
    pa  = oldh['PA']
    cd11 =  0.0006931*sin((pa-0.02)*deg)
    cd12 = -0.0007005*cos((pa+0.02)*deg)
    cd21 =  0.0006931*cos((pa-0.02)*deg)
    cd22 =  0.0007005*sin((pa+0.02)*deg)
    newh=pf.Header()
    newh['SIMPLE']=True
    newh['BITPIX']=8
    newh['NAXIS']=0
    newh['WCSAXES']=2
    newh['CTYPE1']='RA---TAN'
    newh['CTYPE2']='DEC--TAN'
    newh['CRVAL1']=ra
    newh['CRVAL2']=dec
    newh['PA']=pa
    newh['CRPIX1']=600
    newh['CRPIX2']=577
    newh['CUNIT1']='deg     '
    newh['CUNIT2']='deg     '
    newh['CD1_1']=cd11
    newh['CD1_2']=cd12
    newh['CD2_1']=cd21
    newh['CD2_2']=cd22
    newh.add_comment('Simulated astrometry.net WCS based on LE1 FITS keywords')
    newh.tofile(tbd)
  except:
    print(tbd[:-3]+'fits','not present, so not making WCS file.')


