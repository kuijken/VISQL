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

catname=argv[1]

x,y,f,fe,r,a,b,pa,id,flg,rad=np.loadtxt(catname,unpack=True)

if catname[-3:]=='.gz':
    catname=catname[:-3]

bmin=1.3
b13=b>bmin  # filter against cosmic rays, which are thin.
f13=f[b13]
r13=r[b13]

def chimney(r,f,rc=1.1,fc=3e5,rwidth=0.1,fwidth=0.15):
    smokeup=f < (r/rc)**2 * fc *(1+fwidth)
    smokelo=f > (r/rc)**2 * fc /(1+fwidth)
    stackri=r < rc*(1+rwidth)
    stackle=r > rc/(1+rwidth)
    stacklo=f > 1e5 * rc**2 /1.1**2     # modified from 1e5 to give sloped limit
    diaghi =f < (r/rc)**4 * fc
    return (smokeup & smokelo & diaghi & stackle) | (stackri & stackle & stacklo & diaghi)

counthi=0
for rc in np.logspace(0,1,20):
    for fc in np.logspace(0,0.5,8)*3e5*rc**2:     # allow smoke to be higher
        count=chimney(r13,f13,rc,fc).sum()
        if count>counthi:
            counthi=count
            rcbest=rc
            fcbest=fc
#            print (rcbest,fcbest,count)

rmin=rcbest/1.15
rmax=rcbest*1.15
fmax=fcbest/2
fmin=fmax/10
rmed=np.median(r13[(f13>fmin) & (f13<fmax) & (r13<rmax) & (r13>rmin)])
rmin=rmed/1.1
rmax=rmed*1.1
rmed=np.median(r13[(f13>fmin) & (f13<fmax) & (r13<rmax) & (r13>rmin)])
rmin=rmed/1.1
rmax=rmed*1.1

plt.plot(r13,f13,'k,')
ch=chimney(r13,f13,rcbest,fcbest)
plt.plot(r13[ch],f13[ch],'r,')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.5,20)
plt.ylim(1e3,1e7)
plt.plot([rmin,rmax,rmax,rmin,rmin],[fmin,fmin,fmax,fmax,fmin],'g')
plt.title(catname)
plt.xlabel('FLUX RADIUS [pix]')
plt.ylabel('FLUX_AUTO [adu]')
figurename=catname[:-4]+"_chimney2.png"
plt.savefig(figurename)
plt.close()

print('-v fmin=%f -v fmax=%f -v rmin=%f -v rmax=%f -v bmin=%f' % (fmin,fmax,rmin,rmax,bmin))
