from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
from sys import argv
import json
import datetime

rootname=argv[1][:42]
dateobs=datetime.datetime.fromisoformat(
    rootname[2:6]+'-'+rootname[6:8]+'-'+rootname[8:10]+'T'+
    rootname[11:13]+':'+rootname[13:15]+':'+rootname[15:17])

tshutter=0.   # could add shutter delay here

astrom,flux=[x.split() for x in open(rootname+'_x_gaia_dxdy.txt','r')][2:]

mcat,gcat=astrom[:2]

N=int(astrom[2])
sig,dxm,dym,ddxdx,ddydx,ddxdy,ddydy=[float(x) for x in astrom[3:]]

bp1,bp2,bp3,frat1,frat2,frat3,texp=[float(x) for x in flux[2:9]]
texpratio=(89.52+tshutter)/(texp+tshutter)     # ratio of exposure time to fiducial short science

mtab=Table.read(mcat,format='ascii')
gtab=Table.read(gcat,format='ascii')

xm,ym=mtab['XWORLD'],mtab['YWORLD']
xg,yg=gtab['RAV'],gtab['DECV']

# project vis cat to gaia cat using astrometry transformation read in
x=xm+(dxm+(ddxdx*xm+ddxdy*ym)*0.75)/36000
y=ym+(dym+(ddydx*xm+ddydy*ym)*0.75)/36000

# ig is the index in the gaia catalogue of the closest entry to the VIS source
ig=np.array([np.argmin((xx-xg)**2+(yy-yg)**2) for (xx,yy) in zip(x,y)])
disttog=((x-xg[ig])**2+(y-yg[ig])**2)**0.5 * 3600

visok=disttog<2.   # keep all matches better than 2 arcsec
gaiaok=ig[visok]   # gaia index of matching star

visf=mtab['FLUX'][visok]
visx=xm[visok]
visy=ym[visok]
ext=mtab['EXT'][visok]

######## Here do flux correction per quadrant!
led1levelsfile='led1_C_20230802_045303_E_00003003.VIS.bin_02_01_level.txt'
led2levelsfile='led2_C_20230814_082021_C_AC03A900.VIS.bin_04_01_level.txt'
led3levelsfile='led3_C_20230814_085952_C_AC03A980.VIS.bin_01_01_level.txt'
led4levelsfile='led4_C_20230814_085952_C_AC03A980.VIS.bin_03_01_level.txt'
led1ext=np.loadtxt(led1levelsfile,unpack=True,usecols=[1],dtype='str')
led1lev=np.loadtxt(led1levelsfile,unpack=True,usecols=[2])
led2ext=np.loadtxt(led2levelsfile,unpack=True,usecols=[1],dtype='str')
led2lev=np.loadtxt(led2levelsfile,unpack=True,usecols=[2])
led3ext=np.loadtxt(led3levelsfile,unpack=True,usecols=[1],dtype='str')
led3lev=np.loadtxt(led3levelsfile,unpack=True,usecols=[2])
led4ext=np.loadtxt(led4levelsfile,unpack=True,usecols=[1],dtype='str')
led4lev=np.loadtxt(led4levelsfile,unpack=True,usecols=[2])
meanlev1=led1lev.mean()
meanlev2=led2lev.mean()
meanlev3=led3lev.mean()
meanlev4=led4lev.mean()
dsec=led1ext
gain0 = {e:1. for e in dsec}
gain1 = {e:meanlev1/led1lev[led1ext==e][0] for e in dsec}
gain2 = {e:meanlev2/led2lev[led2ext==e][0] for e in dsec}
gain3 = {e:meanlev3/led3lev[led3ext==e][0] for e in dsec}
gain4 = {e:meanlev4/led4lev[led4ext==e][0] for e in dsec}
gainV = json.load(open('EUC_VIS_MDL-GAIN__20231017T184821.468620Z.json','r'))
gainV = {e:gainV[e]/3.4 for e in dsec}

gainchoice={'0':gain0,'1':gain1,'2':gain2,'3':gain3,'4':gain4,'V':gainV}
gainsym={'0':'k','1':'b','2':'c','3':'r','4':'m','V':'g'}

led='3'
gain,gainc=gainchoice[led],gainsym[led]

G=gtab['G'][gaiaok]
bp=gtab['BpRp'][gaiaok]

gflux=10**(-0.4*(G-30))/1.5 # calibrated to short science exposures, roughly
vflux=visf*texpratio*[gain[e] for e in ext]

d1=np.log(frat1)/(bp1-bp2)/(bp1-bp3)
d2=np.log(frat2)/(bp2-bp1)/(bp2-bp3)
d3=np.log(frat3)/(bp3-bp1)/(bp3-bp2)
fratmod= d1*(bp-bp2)*(bp-bp3) + d2*(bp-bp1)*(bp-bp3) + d3*(bp-bp1)*(bp-bp2)
fratmod=np.exp(fratmod)
frat=vflux/gflux/fratmod

#plt.scatter(visx,visy,c=frat,vmin=0.9,vmax=1.1,s=5)
#plt.colorbar()
#plt.show()

# calculate predicted vis flux from fiducial quadratic fit
#  taking into account exposure time

#Collate this info in a new table: ext,visx,visy,G,BpRp,Fvis,frat,Texp

# save as table somehow. pref ascii so can collate and plot in shell scripts easily?

vflux.name='VFLUX'
gflux.name='GFLUX'
frat.name='FRAT'
fratmod.name='FRATMOD'
ttt=Table.Column(0*vflux+texp,name='TEXP')
trat=Table.Column(0*vflux+texpratio,name='TEXPRATIO')

outtab=Table([ext,visx,visy,visf,ttt,G,bp,vflux,gflux,frat,fratmod,trat])
outtab.write(rootname+'_frat.txt',format='ascii',overwrite=True)

#print('%5s %8.5f %8.5f %10.1f  %6.2f  %6.3f %6.3f  %10.1f %10.1f  %8.5f %8.5f %7.4f'
#          % (ext,visx,visy,visf,texp,G,bp,vflux,gflux,frat,fratmod,texpratio))
