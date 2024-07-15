import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import use
import astropy.io.fits as pf

# select proton CRs as having flux between 10000 and 40000 ADU
# plot counts on VIS FP

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')
catfile=argv[1]
x,y,f,fe,rad,a,b,pa,id,flg,rim=np.loadtxt(catfile,unpack=True)

if catfile[-3:]=='.gz':
    catfile=catfile[:-3]
try:
    exptime=pf.open(catfile[:-3]+'fits')[0].header['EXPTIME']
except:
    exptime=90
    print('**** assumed 90 sec exposure time ****')

readtime=72
expfactor= (exptime+readtime/2) / 560   # scale exposure time to 560 sec, for limits of hist2d
    
cosmic=(rad<0.9) & (f<40000) & (f>10000)
nbins=int(0.5+0.8*36000/250)
pixperbin=0.8*36000/nbins
#print(nbins)
aaa=plt.hist2d(x[cosmic],y[cosmic],bins=nbins,range=[[-0.4,0.4],[-0.4,0.4]],
               norm='log',vmin=1*expfactor,vmax=300*expfactor,cmap='jet')
plt.colorbar(label='High-E cosmics / %ix%i pixels' % (pixperbin,pixperbin))
plt.xlabel('VIS X [deg]')
plt.ylabel('VIS Y [deg]')
plt.title(catfile)
plt.savefig(catfile[:-4]+'_protons.png')

counts=aaa[0].flatten()
counts.sort()
top10count=counts[-10]-np.median(counts)
print(catfile[:-4]+'.fits',top10count,exptime+readtime/2,exptime)

