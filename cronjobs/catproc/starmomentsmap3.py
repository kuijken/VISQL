import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from matplotlib import use

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')
plt.rcParams['text.usetex']=False

# Get moments from a catalogue with star measures, eg from starmeaures.py

catname= argv[1]

try:
    extname= np.loadtxt(catname,unpack=True,usecols=0,dtype=str)
    x,y,flx,flxerr,r,rim,xccd,yccd,xc,yc,sumwt,R2,e1,e2,t1,t2,c1,c2 = np.loadtxt(catname,unpack=True,usecols=range(1,19))
except:
    print(catname,'is empty.')
    plt.figure(figsize=(12,9))
    plt.savefig(catname[:-4]+'_map3.png')
    exit()
    
if extname.size<=1:
    print(catname,'is empty.')
    plt.figure(figsize=(12,9))
    plt.savefig(catname[:-4]+'_map3.png')
    exit()

# strip .gz suffix if present
if catname[-3:]=='.gz':
    catname=catname[:-3]
    
exts=np.unique([s[:3] for s in extname])
iext=np.arange(len(exts))
xx=0.*iext
yy=0.*iext
RR2=0.*iext
coma1=0.*iext
coma2=0.*iext
ell1=0.*iext
ell2=0.*iext
tref1=0.*iext
tref2=0.*iext
for i in iext:
    keep=np.char.startswith(extname,exts[i]) & ~np.isnan(R2)
    RR2[i]=np.median(R2[keep])
    coma1[i]=np.median(c1[keep])
    coma2[i]=np.median(c2[keep])
    ell1[i]=np.median(e1[keep])
    ell2[i]=np.median(e2[keep])
    tref1[i]=np.median(t1[keep])
    tref2[i]=np.median(t2[keep])
    xx[i]=np.median(x[keep])
    yy[i]=np.median(y[keep])

plt.figure(figsize=(12,9))

R2min,R2max=np.percentile(RR2,[25,75])
dR2=R2max-R2min
R2min-=dR2
R2max+=dR2

r2d=180./np.pi

coma=(coma1**2+coma2**2)**0.5
cmax=np.percentile(coma,75)
cang=np.arctan2(coma2,coma1)

ell=(ell1**2+ell2**2)**0.5
emax=np.percentile(ell,75)
# halve the azimuth for plotting
eang=np.arctan2(ell2,ell1)/2
ell1,ell2=ell*np.cos(eang),ell*np.sin(eang)

tref=(tref1**2+tref2**2)**0.5
tmax=np.percentile(tref,75)
# divide the azimuth by 3 for plotting
tang=np.arctan2(tref2,tref1)/3
tref1,tref2=tref*np.cos(tang),tref*np.sin(tang)


# set some reasonable ranges for the various aberrations, to make it easier to compare
# plots visually
R2min=min(R2min,1.7)
R2max=max(R2max,2.3)
cmax=max(cmax,0.2)
emax=max(emax,0.07)
tmax=max(tmax,0.4)

plt.subplot(221)
plt.plot(x,y,'k,')
plt.scatter(xx,yy,s=RR2/R2max*50,c=RR2,vmin=R2min,vmax=R2max,cmap='coolwarm')
#plt.xlabel('X')
plt.ylabel('Y')
plt.xlim(-0.4,0.4)
plt.ylim(-0.45,0.45)
plt.colorbar(label='R2')
plt.title(catname[:-3])

plt.subplot(222)
plt.quiver(xx,yy,coma1,coma2,coma,angles='uv',pivot='middle',scale=10*cmax,
               scale_units='width',cmap='Blues',clim=(0,cmax))
plt.quiver(0.36,0.4,cmax,0,cmax,  angles='uv',pivot='middle',scale=10*cmax,
                scale_units='width',cmap='Blues',clim=(0,cmax))
plt.xlim(-0.4,0.4)
plt.ylim(-0.45,0.45)
plt.colorbar(label='Coma')

plt.subplot(223)
plt.quiver(xx,yy,ell1,ell2,ell,angles='uv',pivot='middle',scale=12*emax,
               scale_units='width',
               cmap='Greens',headwidth=1,headlength=0,headaxislength=0,clim=(0,emax))
plt.quiver(0.36,0.4,emax,0,emax,angles='uv',pivot='middle',scale=12*emax,
               scale_units='width',
               cmap='Greens',headwidth=1,headlength=0,headaxislength=0,clim=(0,emax))
plt.xlim(-0.4,0.4)
plt.ylim(-0.45,0.45)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Ellipticity')

plt.subplot(224)
# print the trefoil symbol as 3 rotated quivers
cos,sin=np.cos(2*np.pi/3),np.sin(2*np.pi/3.)
t1scale,t2scale=tmax,0
for rotate in range(3):
    plt.quiver(xx,yy,tref1,tref2,tref,angles='uv',pivot='tail',scale=20*tmax,
                   scale_units='width',cmap='Reds',
                   headwidth=1,headlength=0,headaxislength=0,clim=(0,tmax))
    plt.quiver(0.36,0.4,t1scale,t2scale,tmax,angles='uv',pivot='tail',scale=20*tmax,
                   scale_units='width',cmap='Reds',
                   headwidth=1,headlength=0,headaxislength=0,clim=(0,tmax))
    tref1,tref2=cos*tref1+sin*tref2,-sin*tref1+cos*tref2
    t1scale,t2scale=cos*t1scale+sin*t2scale,-sin*t1scale+cos*t2scale
plt.xlim(-0.4,0.4)
plt.ylim(-0.45,0.45)
plt.xlabel('X')
#plt.ylabel('Y')
plt.colorbar(label='Trefoil')

plt.tight_layout()
plt.savefig(catname[:-4]+'_map3.png')

