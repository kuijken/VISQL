import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import use

use('agg')

catname=argv[1]

x,y,f,fe,r,a,b,pa,id,flg,rad=np.loadtxt(catname,unpack=True)

if catname[-3:]=='.gz':
    catname=catname[:-3]

bmin=1.3
b13=b>bmin  # filter against cosmic rays, which are thin.
f13=f[b13]
r13=r[b13]

def chimney(r,f):
    return (r>1) & (r<2) & (f<3e5) & (f>1e4)


rmin=1
rmax=2
fmax=3e5
fmin=1e4
rmed=1.5

plt.plot(r13,f13,'k,')
ch=chimney(r13,f13)
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
