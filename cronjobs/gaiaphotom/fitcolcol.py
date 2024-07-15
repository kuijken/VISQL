from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
from sys import argv
import json
import datetime
from sys import argv

if len(argv)<2:
    print('Usage: p3 fitcolcol.py <rootname of frat.txt file> exptime [mode]   (default mode is 1)')
    print('  mode=1: filter on 3x3 blocks only')
    print('  mode=2: filter on 3x3 blocks and ccds')
    print('  mode=3: filter on 3x3 blocks, ccds and quadrants')
    exit()
    
root=argv[1][:42]
texp=float(argv[2])
tshutter=0
texpratio=(89.52+tshutter)/(texp+tshutter)  # ratio wrt short science
verbose=False
filtermode=1
if len(argv)>3:
    filtermode=int(argv[3])

t=Table.read(root+'_frat.txt',format='ascii')
if len(t)==0:
    print('No data')
    exit()

# gain tables to correct default table 3 to table V
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

oldgain=np.array([gain3[e] for e in t['EXT']])
newgain=np.array([gainV[e] for e in t['EXT']])
                     
# correct VIS fluxes to gain table V from gain table 3.
t['VFLUX'] *= newgain/oldgain

keep=t['TEXP']==texp

xw=t['XWORLD'][keep]
yw=t['YWORLD'][keep]
ext=t['EXT'][keep]

ccd=np.array([xt[:3] for xt in ext])
ccds=np.unique(ccd)

#filters={
#    'lole':(t['XWORLD']<-0.11) & (t['YWORLD']<-0.13),
#    'upri':(t['XWORLD']> 0.11) & (t['YWORLD']> 0.13)
#    }

filters={
    '1212': (xw<-0.11)             & (yw<-0.13),
    '1234': (xw<-0.11)             & (yw>-0.13) & (yw<0.13),
    '1256': (xw<-0.11)             & (yw> 0.13),
    '3412': (xw>-0.11) & (xw<0.11) & (yw<-0.13),
    '3434': (xw>-0.11) & (xw<0.11) & (yw>-0.13) & (yw<0.13),
    '3456': (xw>-0.11) & (xw<0.11) & (yw> 0.13),
    '5612': (xw> 0.11)             & (yw<-0.13),
    '5634': (xw> 0.11)             & (yw>-0.13) & (yw<0.13),
    '5656': (xw> 0.11)             & (yw> 0.13)
    }

if filtermode>1:
    for c in ccds:
        filters['CCD'+c[0]+c[2]]=(ccd==c)
if filtermode>2:
    for xt in np.unique(ext):
        filters[xt]=(ext==xt)

for filt in filters:
 br=t['BpRp'][keep][filters[filt]]
 frat=np.log(t['VFLUX']/t['GFLUX'])[keep][filters[filt]]

 # fit log(fluxratio) vs bprp (quadratic poly)
 one=np.ones(len(br))
 vec=np.vstack([one,br,br**2])
 coefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([frat]).T,rcond=None)
 rms=(residsq/len(br))**0.5
 fit=coefs[0]+br*(coefs[1]+br*coefs[2])
 ok=(fit-frat)**2<9*rms**2
 #print('clipping',(ok==False).sum(),'Gaia Bp-Rp stars')
 # Iterate out outliers
 brfit=br[ok]
 one=np.ones(len(brfit))
 vec=np.vstack([one,brfit,brfit**2])
 coefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([frat[ok]]).T,rcond=None)
 meanbr=brfit.mean()

 xx=np.linspace(0.5,3,26)
 yy=coefs[0]+xx*(coefs[1]+xx*coefs[2])
 #save fitted value at bp-rp=0.75,1.75,2.75
 bpref=np.linspace(0.75,2.75,3)
 fratref=coefs[0]+bpref*(coefs[1]+bpref*coefs[2])

 print(root+'_stars2_mom.cat','%sVIS/Gaia 0.75 1.75 2.75 %7.4f %7.4f %7.4f %9.2f' %
          ((filt,)+tuple(np.exp(fratref))+(texp,)))

 
