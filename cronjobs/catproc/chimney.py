import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import astropy.io.fits as pf
from matplotlib import use

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')

bamin=0.5
rfieldmin=0.1
for cat in argv[1:]:
    print('Plotting data from',cat)
    x,y,f,fe,r,a,b,pa,id,flg,fieldr=np.loadtxt(cat,unpack=True)
    exptime=pf.open(cat[:-6]+'fits')[0].header['EXPTIME']
    sub=(r>0.5) & (fieldr < rfieldmin) & (b>bamin*a)
    plt.figure(figsize=(8,6))
    plt.scatter(r[sub],25.67-2.5*np.log10(f[sub]*3.4/exptime),s=(b/a)[sub],c=1-(b/a)[sub],cmap='jet',vmax=1-bamin,vmin=0)
    plt.colorbar(label='ellipticity 1-b/a')
    plt.xlim(0,5)
    plt.ylim(26,16)
    plt.xlabel('FLUXRADIUS [pix]')
    plt.ylabel('MAG_AUTO [approx]')
    plt.title(cat[:-7] + ' - inner %.2f deg' % rfieldmin)
    plt.grid()
    plt.savefig(cat[:-7]+'_chimney4.png')
    plt.close()
    

