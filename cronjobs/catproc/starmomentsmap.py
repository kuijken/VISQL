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

catname= argv[1]

try:
    extname= np.loadtxt(catname,unpack=True,usecols=0,dtype=str)
    x,y,flx,flxerr,r,rim,xccd,yccd,xc,yc,sumwt,R2,e1,e2,t1,t2,c1,c2 = np.loadtxt(catname,unpack=True,usecols=range(1,19))
except:
    print(catname,'is empty.')
    exit()

if len(extname)==0:
    print(catname,'is empty.')
    exit()
       
if catname[-3:]=='.gz':
    catname=catname[:-3]

plt.figure(figsize=(10,10))

# subsample the data if there are too many points to plot.
if len(x)>1000:
    keep=np.random.rand(len(x)) < 1000/len(x)
    x=x[keep]
    y=y[keep]
    c1=c1[keep]
    c2=c2[keep]
    e1=e1[keep]
    e2=e2[keep]
    t1=t1[keep]
    t2=t2[keep]
    R2=R2[keep]
    rim=rim[keep]

cmax=max(np.percentile(np.abs(c1),75),np.percentile(np.abs(c2),75))
if np.isnan(cmax):    cmax=1
cmax=0.15
plt.subplot(441)
plt.plot(c1,c2,'b.',markersize=1)
plt.xlim(-2*cmax,2*cmax)
plt.ylim(-2*cmax,2*cmax)
plt.grid()
plt.xlabel('Coma1')
plt.ylabel('Coma2')
plt.subplot(442)
plt.plot(x,c1,'b.',markersize=1)
plt.plot(x,c2,'r.',markersize=1)
plt.grid()
plt.xlabel('X')
plt.xlim(-0.35,0.35)
plt.ylim(-2*cmax,2*cmax)
plt.title(catname)
plt.subplot(443)
plt.plot(y,c1,'b.',markersize=1)
plt.plot(y,c2,'r.',markersize=1)
plt.xlabel('Y')
plt.xlim(-0.4,0.4)
plt.ylim(-2*cmax,2*cmax)
plt.subplot(444)
plt.plot(rim,c1,'b.',markersize=1)
plt.plot(rim,c2,'r.',markersize=1)
plt.xlabel('R')
plt.xlim(0,0.55)
plt.ylim(-2*cmax,2*cmax)

emax=max(np.percentile(np.abs(e1),75),np.percentile(np.abs(e2),75))
if np.isnan(emax):    emax=1
emax=0.05
plt.subplot(445)
plt.plot(e1,e2,'b.',markersize=1)
plt.xlim(-2*emax,2*emax)
plt.ylim(-2*emax,2*emax)
plt.grid()
plt.xlabel('ell1')
plt.ylabel('ell2')
plt.subplot(446)
plt.plot(x,e1,'b.',markersize=1)
plt.plot(x,e2,'r.',markersize=1)
plt.xlabel('X')
plt.xlim(-0.35,0.35)
plt.ylim(-2*emax,2*emax)
plt.subplot(447)
plt.plot(y,e1,'b.',markersize=1)
plt.plot(y,e2,'r.',markersize=1)
plt.xlabel('Y')
plt.xlim(-0.4,0.4)
plt.ylim(-2*emax,2*emax)
plt.subplot(448)
plt.plot(rim,e1,'b.',markersize=1)
plt.plot(rim,e2,'r.',markersize=1)
plt.xlabel('R')
plt.xlim(0,0.55)
plt.ylim(-2*emax,2*emax)

tmax=max(np.percentile(np.abs(t1),75),np.percentile(np.abs(t2),75))
if np.isnan(tmax):    tmax=1
tmax=0.25
plt.subplot(449)
plt.plot(t1,t2,'b.',markersize=1)
plt.xlim(-2*tmax,2*tmax)
plt.ylim(-2*tmax,2*tmax)
plt.grid()
plt.xlabel('Tref1')
plt.ylabel('Tref2')
plt.subplot(4,4,10)
plt.plot(x,t1,'b.',markersize=1)
plt.plot(x,t2,'r.',markersize=1)
plt.xlabel('X')
plt.xlim(-0.35,0.35)
plt.ylim(-2*tmax,2*tmax)
plt.subplot(4,4,11)
plt.plot(y,t1,'b.',markersize=1)
plt.plot(y,t2,'r.',markersize=1)
plt.xlabel('Y')
plt.xlim(-0.4,0.4)
plt.ylim(-2*tmax,2*tmax)
plt.subplot(4,4,12)
plt.plot(rim,t1,'b.',markersize=1)
plt.plot(rim,t2,'r.',markersize=1)
plt.xlabel('R')
plt.xlim(0,0.55)
plt.ylim(-2*tmax,2*tmax)

R2min,R2max=np.percentile(R2,[25,75])
dR2=R2max-R2min
R2min-=dR2
R2max+=dR2
if np.isnan(dR2):  R2min,R2max=1,10
R2min,R2max=1.7,2.4
#plt.subplot(4,4,13)
#plt.plot(x,y,'k.',markersize=1)
#plt.xlabel('X')
#plt.ylabel('Y')
plt.subplot(4,4,14)
plt.plot(x,R2,'k.',markersize=1)
plt.xlabel('X')
plt.ylabel('R2')
plt.xlim(-0.35,0.35)
plt.ylim(R2min,R2max)
plt.subplot(4,4,15)
plt.plot(y,R2,'k.',markersize=1)
plt.xlabel('Y')
plt.xlim(-0.4,0.4)
plt.ylim(R2min,R2max)
plt.subplot(4,4,16)
plt.plot(rim,R2,'k.',markersize=1)
plt.xlabel('R')
plt.xlim(0,0.55)
plt.ylim(R2min,R2max)

plt.tight_layout()
plt.savefig(catname[:-4]+'.png')
