import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import os
from matplotlib import use

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')

for catfile in argv[1:]:
    if catfile[-3:]=='.gz':
        pngfile = catfile[:-7]+'_rf.png'
    else:
        pngfile = catfile[:-4]+'_rf.png'            
    if os.path.exists(pngfile):
        #print(pngfile,'already exists.')
        pass
    else:
        f,r,flg=np.loadtxt(catfile, usecols=(2,4,9), unpack=True)
        plt.hist2d(r,np.log10(f),bins=100,range=[[0,3],[2.5,5.5]],cmap='jet',norm='log',vmin=1e2,vmax=1e4)
        plt.colorbar()
        plt.xlabel('Flux Radius [pix]')
        plt.ylabel('log10(flux) [ADU]')
        plt.title(catfile)
        if catfile[-3:]=='.gz':
            catfile=catfile[:-3]
        plt.savefig(pngfile)
        plt.close()
        print(pngfile,'created.')

# density up to 50k/pix in 289s exposure 
