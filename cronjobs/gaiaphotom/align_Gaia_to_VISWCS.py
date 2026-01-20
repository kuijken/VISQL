import astropy.io.fits as pf
import astropy.io.votable as vot
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
from sys import argv
from os import environ
import match


plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')

tiledir=environ['viskom']+'/SCIENCE/'

catname=argv[1]
rac,decc,pa=[float(x) for x in argv[2:5]]

#tilename='EUC_VIS_SWL-STK-000-000000-0000000__20231121T232926.310495Z.cat'
#tile=pf.open(tiledir+tilename)
tile=np.loadtxt(tiledir+catname,unpack=True)
if not len(tile):
    exit()
    
tilecat={}
tilecat['RAVIS']=tile[0]
tilecat['DECVIS']=tile[1]
tilecat['VISFLUX']=tile[2]
tilecat['VISFLUXERR']=tile[3]
tilecat['FLUX_RADIUS']=tile[4]
tilecat['A_IMAGE']=tile[5]
tilecat['B_IMAGE']=tile[6]
tilecat['THETA_IMAGE']=tile[7]
tilecat['NUMBER']=tile[8]
tilecat['FLAGS']=tile[9]
del tile

if catname[-3:]=='.gz': catname=catname[:-3]

crac =np.cos(np.pi/180*rac)
srac =np.sin(np.pi/180*rac)
cdecc=np.cos(np.pi/180*decc)
sdecc=np.sin(np.pi/180*decc)
cpa  =np.cos(np.pi/180*(pa+90))
spa  =np.sin(np.pi/180*(pa+90))

ratile= tilecat['RAVIS']
dectile=tilecat['DECVIS']

####now read in Gaia patch, cut down to tile.
# gaia=pf.open('GaiaDR3xSDSSDR16SelfCal.fits')[1].data
# visragaia =gaia['ra']
# visdecgaia=gaia['dec']
# G=gaia['phot_g_mean_mag']
# Gerr=1.086/gaia['phot_g_mean_flux_over_error']    # (1.086=2.5/ln(10))
# bprp=gaia['bp_rp']
# idgaia=gaia['source_id']
###

### or query from vizier
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
Vizier.ROW_LIMIT=-1
Vizier.catalog='I/355/gaiadr3'
try:
    gaia=Vizier.query_region(coord.SkyCoord(ra=rac,dec=decc,unit=(u.deg,u.deg)),radius=0.6*u.deg)[0]
except:
    print('Vizier query returned no results for',catname,'around RA,DEC=',rac,decc)
    exit()

visragaia =gaia['RA_ICRS']
visdecgaia=gaia['DE_ICRS']
G=gaia['Gmag']
Gerr=1.086*gaia['e_FG']/gaia['FG']    # (1.086=2.5/ln(10))
bprp=gaia['BP-RP']
idgaia=gaia['Source']
###

# rotate Gaia into coordinates of the tile's VIS WCS
cra= np.cos(np.pi/180*(visragaia-rac))
sra= np.sin(np.pi/180*(visragaia-rac))
cdec=np.cos(np.pi/180*visdecgaia)
sdec=np.sin(np.pi/180*visdecgaia)
# rotate by racc, then -decc, then PA
x=cra*cdec*cdecc + sdec*sdecc
y=sra*cdec*cpa+cra*cdec*sdecc*spa-sdec*cdecc*spa
z=sra*cdec*spa-cra*cdec*sdecc*cpa+sdec*cdecc*cpa
# VIS WCS has flipped RA (which can be plotted L to R therefore)
visragaia= -np.arctan2(y,x) * 180/np.pi
visdecgaia= np.arcsin(z)    * 180/np.pi

# fine tuning step to take out higher-order distortions approximately
visragaia   += (visragaia-0.05)**2*550/36000 +(visdecgaia+0.05)**2*200/36000  # quadratic RA dependence of residual ~ 100 pix amplitude
visdecgaia  += -(visdecgaia*450/36000) + (visdecgaia*visragaia)*350/36000

ontile=(-0.7<visragaia) & (0.7>visragaia) & (-0.7<visdecgaia) & (0.7>visdecgaia)
ontile &= (Gerr<0.02)
ontile &= G<21
ontile &= ~bprp.mask
visragaia =visragaia[ontile]
visdecgaia=visdecgaia[ontile]
fig,ax=plt.subplots()
cc=ax.scatter(visragaia,visdecgaia,
                  s=np.maximum(0,7*(15-G[ontile])),c=G[ontile],
                  vmin=15,vmax=8,cmap='gray')
fig.colorbar(cc,label='Gaia G magnitude')
ax.set_xlim(-0.5,0.5)
ax.set_ylim(-0.5,0.5)
#ax.plot([-0.35,0.35,0.35,-0.35,-0.35],[-0.38,-0.38,0.38,0.38,-0.38],'r')
w=(0.1151-0.0017)
for iccd in range(6):
    xll=iccd*(0.1186-0.0017)-.3484
    for jccd in range(6):
        yll=jccd*(0.1405-0.0087)-.3866
        ax.plot([xll,xll+w,xll+w,xll,xll],[yll,yll,yll+w,yll+w,yll],'r',linewidth=0.5)
        ax.plot([xll,xll+w],[yll+w/2,yll+w/2],'r:',linewidth=0.5)
        ax.plot([xll+w/2,xll+w/2],[yll,yll+w],'r:',linewidth=0.5)
        
arrowl=0.39/max(np.abs(spa),np.abs(cpa))
ax.arrow(0,0,spa*arrowl,cpa*arrowl,width=0.003,head_width=0.03)
ax.set_title('Gaia stars in '+catname)
ax.set_xlabel('RA [VIS instr. WCS]')
ax.set_ylabel('DEC [VIS instr. WCS]')
ax.set_aspect('equal')
plt.text(0,0.45,'RA %7.3f DEC %7.3f PA %6.2f' % (rac,decc,pa),bbox=dict(facecolor='red', alpha=0.1),ha='center')
plt.savefig(catname[:-4]+'_gaia.png')
plt.clf()

G=G[ontile]
Gerr=Gerr[ontile]
bprp=bprp[ontile]
idgaia=idgaia[ontile]

one=np.ones(len(visragaia))

q=(visragaia**2<0.4**2) & (visdecgaia**2<0.4**2)
one=one[q]

# addition: include also the original Gaia RA and DEC in the output _gaia.cat file
#gtab=Table([visragaia[q],visdecgaia[q],G[q],Gerr[q],one,one,one,one,idgaia[q],one,bprp[q],gaia['RA_ICRS'][ontile][q],gaia['DE_ICRS'][ontile][q]],
#               names=['RAV','DECV','G','Gerr','R','A','B','PA','ID','FLG','BpRp','RA_ICRS','DE_ICRS'])
gtab=Table([visragaia[q],visdecgaia[q],G[q],Gerr[q],one,one,one,one,idgaia[q],one,bprp[q]],
               names=['RAV','DECV','G','Gerr','R','A','B','PA','ID','FLG','BpRp'])
gtab.write(catname[:-4]+'_gaia.cat',format='ascii.commented_header',overwrite=True)


exit()

#then run matchpos to find the matches
matchid,sep=match.matchpos(ratile,dectile,ragaia,decgaia,1./3600)
   # matchid are the indices in tile catalogue of matched stars; -1 means no match within 1 arcsec

#get all photometry for these stars. For Gaia keep those with matchid>-1
ragaia=ragaia[matchid>-1]
decgaia=decgaia[matchid>-1]
G=G[ontile][matchid>-1]
Gerr=Gerr[ontile][matchid>-1]
g=gaia['gmag'][ontile][matchid>-1]
r=gaia['rmag'][ontile][matchid>-1]
i=gaia['imag'][ontile][matchid>-1]
idgaia=idgaia[ontile][matchid>-1]
# and for tile stars just get the indices that are nonnegative
matchid=matchid[matchid>-1]
visflux=tilecat['VISFLUX'][matchid]
vis=-2.5*np.log10(visflux)
rad=tilecat['FLUX_RADIUS'][matchid]

gi=g-i
ri=r-i
Vr=vis-r
Vi=vis-i
# now attempt to locate the SLR magnitudes.
# first locate the centre of the g-i distribution - this is roughly the knee in the r-G v g-i diagram

plt.scatter(gi,(0.3*Vr+0.7*Vi),s=1)
#plt.colorbar(label='FLUX_RADIUS')
plt.xlabel('g-i (SDSS)')
plt.ylabel('VIS-r(SDSS)')
plt.xlim(0.3,2)
plt.ylim(-24.,-23.5)
plt.title(catname)
plt.savefig(catname[:-4]+'.Gaia.png')



#to query gaia on the fly, use sth like (but with correct catalogue name)
#from astroquery.vizier import Vizier
#import astropy.coordinates as coord
#import astropy.units as u
#Vizier.ROW_LIMIT=-1
#Vizier.catalog='I/355/gaiadr3'
#gaia=Vizier.query_region(coord.SkyCoord(ra=268,dec=65,unit=(u.deg,u.deg)),radius=1*u.deg)[0]
#visragaia =gaia['RA_ICRS']
#visdecgaia=gaia['DE_ICRS']
#G=gaia['Gmag']
#Gerr=1.086*gaia['FG']/gaia['e_FG']    # (1.086=2.5/ln(10))
#bprp=gaia['BP-RP']
#idgaia=gaia['Source']
