import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import use
import astropy.io.fits as pf

# select X rays as having flux between 400 and 600 ADU
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
    
cosmic=(rad<0.9) & (f<1200) & (f>400)
nbins=int(0.5+0.8*36000/250)
pixperbin=0.8*36000/nbins
#print(nbins)
aaa=plt.hist2d(x[cosmic],y[cosmic],bins=nbins,range=[[-0.4,0.4],[-0.4,0.4]],
               norm='log',vmin=1*expfactor,vmax=300*expfactor,cmap='jet')
plt.colorbar(label='Cosmics / %ix%i pixels' % (pixperbin,pixperbin))
plt.xlabel('VIS X [deg]')
plt.ylabel('VIS Y [deg]')
plt.title(catfile)
plt.savefig(catfile[:-4]+'_cosmiX.png')
plt.clf()

counts=aaa[0].flatten()
counts.sort()
top10count=counts[-10]-np.median(counts)
print(catfile[:-4]+'.fits',top10count,exptime+readtime/2,exptime)

# make cosmics plot while here.

cosmic=(rad<0.9) & (f<3000)
nbins=int(0.5+0.8*36000/200)
pixperbin=0.8*36000/nbins
#print(nbins)
plt.hist2d(x[cosmic],y[cosmic],bins=nbins,range=[[-0.4,0.4],[-0.4,0.4]],
               norm='log',vmin=10,vmax=2000,cmap='jet')
plt.colorbar(label='Cosmics / %ix%i pixels' % (pixperbin,pixperbin))
plt.xlabel('VIS X [deg]')
plt.ylabel('VIS Y [deg]')
plt.title(catfile)
plt.savefig(catfile[:-4]+'_cosmics.png')
plt.clf()

# make protons plot while here.

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
plt.clf()


counts=aaa[0].flatten()
counts.sort()
top10count=counts[-10]-np.median(counts)
#print(catfile[:-4]+'.fits',top10count,exptime+readtime/2,exptime)

# make rfdiagram while here.

pngfile=catfile[:-4]+'_rf.png'
plt.hist2d(rad,np.log10(f),bins=100,range=[[0,3],[2.5,5.5]],cmap='jet',norm='log',vmin=1e2,vmax=1e4)
plt.colorbar()
plt.xlabel('Flux Radius [pix]')
plt.ylabel('log10(flux) [ADU]')
plt.title(catfile)
plt.savefig(pngfile)
plt.clf()
print(pngfile,'created.')

plt.close()
