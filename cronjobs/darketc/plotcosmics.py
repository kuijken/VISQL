import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import use

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

cosmic=(rad<0.9) & (f<3000)
nbins=int(0.5+0.8*36000/200)
#print(nbins)
plt.hist2d(x[cosmic],y[cosmic],bins=nbins,range=[[-0.4,0.4],[-0.4,0.4]],
               norm='log',vmin=10,vmax=2000,cmap='jet')
plt.colorbar(label='Cosmics / 250x250 pixels')
plt.xlabel('VIS X [deg]')
plt.ylabel('VIS Y [deg]')
plt.title(catfile)
plt.savefig(catfile[:-4]+'_cosmics.png')

