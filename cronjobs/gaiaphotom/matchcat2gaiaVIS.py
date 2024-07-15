import matplotlib.pyplot as plt
import numpy as np
from os.path import basename
from sys import argv
from astropy.table import Table
import json

from matplotlib import use
plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'
use('Agg')

def fitfluxslope(f1,f2):
    '''
    Find best-fit slope f2=S.f1 for flux measurements assuming photon noise only
    minimise sum[  (f1_i-f_i)^2/f1_i + (f2_i- S.f_i)^2/f2_i ] over all fluxes f_i as well as the slope parameter S.
    First minimise over individual fluxes f_i, for a given S:
    f = f1 f2 (1+S) /(f2 + S^2 f1)   for each i
    Then minimize over S, substituting these f_i:
    Min. sum[ ... or just calculate average of log(f1/f2) with proper error analysis?
    '''



if len(argv)<3:
    print('Usage: p3 matchmomcat2gaia.py xx_gaia.cat [exptime [verbose]')
    print('        xx_gaia.cat comes from align_Gaia_to_VISWCS.py')
    exit()

try:
    momcatvis=argv[1]
except:
    exit()

gaiavis=momcatvis[:-8]+'_gaia.cat'

tshutter=0   # could include this if you want - eg 1.5s average shutter delay
texp=89.52  # short science
texpratio=1
rotate=0
verbose=False
if len(argv)>2:
    texp=float(argv[2])
    texpratio=(89.52+tshutter)/(texp+tshutter)
if len(argv)>3:
    verbose=True

mcat=Table.read(momcatvis,format='ascii')
gcat=Table.read(gaiavis,format='ascii')

######## Here do flux correction per quadrant!
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
# pick an LED to derive the gains from, mark colour to distinguish
#gain,gainc=gain0,'k'  
#gain,gainc=gain1,'b'  
#gain,gainc=gain3,'r'
#gain,gainc=gain4,'m'
#gain,gainc=gainV,'g'

led='3'
gain,gainc=gainchoice[led],gainsym[led]

########


x1,y1,f1=mcat['XWORLD'],mcat['YWORLD'],mcat['FLUX'] * texpratio    # apply flux correction factor
ext=mcat['EXT']
print(f1[:10])
f1 *= [gain[e] for e in ext]
print(f1[:10])
x2,y2,f2=gcat['RAV'],gcat['DECV'],10**(-0.4*(gcat['G']-30))/1.5    # zpt is calibrated to short exposures
bprp=gcat['BpRp']
R2=mcat['R2']

c,s=np.cos(np.pi/180*rotate),np.sin(np.pi/180*rotate)
x2,y2=x2*c-y2*s,x2*s+y2*c

xrat=20 #recommend 1.1 for vis x vis catalogues
n1=len(f1)
n2=len(f2)
print(momcatvis,": %i stars found" %n1) 
print(gaiavis,": %i stars found" %n2) 
# only keep objects in the largest catalogue that are within a factor of xrat of
#   the flux range in the other cat, allowing for the exp.time ratio texpratio
if (n1>n2):
    f2min=f2.min()
    f2max=f2.max()
    keep1=(f1>f2min/xrat) & (f1<f2max*xrat)
    x1=x1[keep1]
    y1=y1[keep1]
    f1=f1[keep1]
    R2=R2[keep1]
else:
    f1min=f1.min()
    f1max=f1.max()
    keep2=(f2>f1min/xrat) & (f2<f1max*xrat)
    x2=x2[keep2]
    y2=y2[keep2]
    f2=f2[keep2]
    bprp=bprp[keep2]
n1=len(f1)
n2=len(f2)
print(momcatvis,": %i stars kept" %n1) 
print(gaiavis,": %i stars kept" %n2)
while np.maximum(n1,n2)>5000:
    f110,f190=np.percentile(f1,[10,90])
    f210,f290=np.percentile(f2,[10,90])
    f1max=np.maximum(f190,f290)
    f1min=np.minimum(f110,f210)
    keep1=(f1<f1max) & (f1>f1min)
    keep2=(f2<f1max) & (f2>f1min)
    x1=x1[keep1]
    y1=y1[keep1]
    f1=f1[keep1]
    x2=x2[keep2]
    y2=y2[keep2]
    f2=f2[keep2]
    bprp=bprp[keep2]
    R2=R2[keep1]
    n1=len(f1)
    n2=len(f2)
    print(momcatvis,": %i stars kept" %n1,f1min,f1max) 
    print(gaiavis,": %i stars kept" %n2,f1min,f1max)

dx=np.array([xx2-xx1 for xx2 in x2 for xx1 in x1])
dy=np.array([yy2-yy1 for yy2 in y2 for yy1 in y1])

dtol=0.007 # matching tolerance [deg] (recomment 0.001 for matchng 2 vis cats)
win,winzoom=0.4,5
#win,winzoom=0.01,3
keep=(dx>-win) & (dx<win) & (dy>-win) & (dy<win)
dx=dx[keep]
dy=dy[keep]
for iter in range(2):
 h2=np.histogram2d(dx,dy,bins=200)
 pk=np.unravel_index(h2[0].argmax(),h2[0].shape)
 xoff,yoff=h2[1][pk[0]],h2[2][pk[1]]
 if verbose:
    print(win,xoff,yoff)
 win/=winzoom
 keep=(dx>xoff-win) & (dx<xoff+win) & (dy>yoff-win) & (dy<yoff+win)
 dx=dx[keep]
 dy=dy[keep]
 if verbose:
    plt.imshow(h2[0].T)
    plt.show()
if verbose:
    plt.close()

i2=np.zeros(len(x1),dtype=int)

for i1 in range(len(x1)):
    xx=x1[i1]
    yy=y1[i1]
    d2=(x2-xoff-xx)**2 + (y2-yoff-yy)**2
    ibest=np.argmin(d2)
    if d2[ibest]<dtol**2:
        i2[i1]=ibest
    else:
        i2[i1]=-1

ii2=i2[i2>=0]
x1match=x1[i2>=0]
y1match=y1[i2>=0]
f1match=f1[i2>=0]
x2match=x2[ii2]
y2match=y2[ii2]
f2match=f2[ii2]
bprp=bprp[ii2]
R2=R2[i2>=0]

frat=f2match/f1match

xm=(x1match+x2match)/2
ym=(y1match+y2match)/2
one=np.ones(len(xm))
keep=np.ones(len(one),dtype=bool)
vec=np.vstack([one,xm,ym])
dx=(x2match-x1match)
dy=(y2match-y1match)
dxy=np.vstack([dx,dy])
if verbose:
    print('iter keep  sigxy      dx        dy     ddx/dx   ddy/dx   ddx/dy   ddy/dy')
for iter in range(25):
    coef,residsq,rank,s=np.linalg.lstsq((vec*keep).T,(dxy*keep).T,rcond=None)
    sigxy=(residsq.sum()/keep.sum()/2)**0.5
    keep=((np.dot(vec.T,coef)-dxy.T)**2).sum(axis=1) < 2*(2.5*sigxy)**2
    if verbose:
        print('%2i %5i %7.3f %9.3f %9.3f %8.2f %8.2f %8.2f %8.2f' % (
            (iter,keep.sum(),sigxy*36000)+
            tuple(coef[0,:]*36000)+
            tuple(coef[1:,:].flatten()*36000*0.75) ) )

fwt=(1/f1match+1/f2match)**(-1) # (1/sigma[ratio]**2 = (1/f1 + 1/f2) if Poisson)
logmeanfluxratio=(np.log(f2match[keep]/f1match[keep])*fwt[keep]).sum()/fwt[keep].sum()
meanfluxratio=np.exp(logmeanfluxratio)

dx0pix=coef[0,0]*36000
dy0pix=coef[0,1]*36000
delplot=(dx[keep].std()**2+dy[keep].std()**2)**0.5*4 * 36000

plt.figure(figsize=(12,7))

plt.subplot(231)
plt.plot(x1match,36000*dx,'k.')
plt.plot(x1match[keep],36000*dx[keep],'r.')
#plt.xlabel('x1 [deg]')
plt.ylabel('x2-x1 [pix]')
plt.title(basename(momcatvis))
plt.ylim(dx0pix-delplot,dx0pix+delplot)
plt.xlim(-0.41,0.41)

plt.subplot(234)
plt.plot(x1match,36000*dy,'k.')
plt.plot(x1match[keep],36000*dy[keep],'r.')
plt.xlabel('x1 [deg]')
plt.ylabel('y2-y1 [pix]')
plt.ylim(dy0pix-delplot,dy0pix+delplot)
plt.xlim(-0.41,0.41)

plt.subplot(232)
plt.plot(y1match,36000*dx,'k.')
plt.plot(y1match[keep],36000*dx[keep],'r.')
#plt.xlabel('y1 [deg]')
#plt.ylabel('x2-x1 [pix]')
plt.ylim(dx0pix-delplot,dx0pix+delplot)
if (rotate != 0):
    plt.title('Rotation 2 to 1: %6.1f deg'%rotate)
plt.xlim(-0.41,0.41)

plt.subplot(235)
plt.plot(y1match,36000*dy,'k.')
plt.plot(y1match[keep],36000*dy[keep],'r.')
plt.xlabel('y1 [deg]')
#plt.ylabel('y2-y1 [pix]')
plt.ylim(dy0pix-delplot,dy0pix+delplot)
plt.xlim(-0.41,0.41)

ddx=dx-coef[0,0]
ddy=dy-coef[0,1]
plt.subplot(233)
plt.quiver(xm[keep],ym[keep],ddx[keep],ddy[keep],(ddx**2+ddy**2)[keep]**0.5*36000,cmap='jet')
plt.xlim(-0.4,0.4)
plt.ylim(-0.4,0.4)
plt.colorbar(label='shift rel. to mean [pix]')
plt.title(basename(gaiavis))

# plt.subplot(233)
# plt.scatter(x1match,y1match,c=frat,vmin=0.9,vmax=1.1,cmap='jet')
# plt.colorbar(label='Flux2/Flux1')
# plt.title(cat2)

plt.subplot(236)
#plt.plot(f2match/1e4,f1match/1e4,'k.')
#plt.plot(f2match[keep]/1e4,f1match[keep]/1e4,'r.')
plt.scatter(f2match[keep]/1e4,f1match[keep]/1e4,c=bprp[keep],vmin=0,vmax=3,linewidths=0,cmap='copper')
plt.colorbar(label='Bp-Rp')
plt.plot([f2.min()/1e4,f2.max()/1e4],[f2.min()/1e4,f2.max()/1e4],'r')
plt.grid()
plt.xlabel('GAIA FLUX (arb. units)')
plt.ylabel('VIS FLUX (1E4 ADU)')
# optimal plotting range with origin included
#plt.xlim(0,x1)
#plt.ylim(0,x1*meanfluxratio)
# force limits to correspond to gaia range and scaled VIS
plt.xlim(0,8)
plt.ylim(0,8)
x0,x1=plt.xlim()
plt.plot([0,x1],[0,x1/meanfluxratio],'k',label='slope %.3e' % (1/meanfluxratio))
plt.legend()

#plt.tight_layout()

momcatvis=basename(momcatvis)
gaiavis=basename(gaiavis)

plotname=momcatvis[:-8] + "_x_gaia_dxdy.png"

plt.savefig(plotname)
plt.clf()

lf2=30-2.5*np.log10(1.5*f2match)
keepbr=keep & ~ np.isnan(bprp) & (lf2<19.5)    # fainter stars have dubious BpRp?
br=bprp[keepbr]
R2=R2[keepbr]
frat=np.log(f1match[keepbr]/f2match[keepbr])
lf2=lf2[keepbr]
xm=xm[keepbr]
ym=ym[keepbr]

# fit log(fluxratio) vs bprp (quadratic poly)
one=np.ones(len(br))
vec=np.vstack([one,br,br**2])
coefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([frat]).T,rcond=None)
rms=(residsq/len(br))**0.5
fit=coefs[0]+br*(coefs[1]+br*coefs[2])
ok=(fit-frat)**2<9*rms**2
print('clipping',(ok==False).sum(),'Gaia Bp-Rp stars')
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

# fit R2 vs xm,ym (quadratic) and bprp (linear)
vec=np.vstack([np.ones(len(xm)),xm,ym,xm**2,xm*ym,ym**2,br])
#xycoefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([R2]).T,rcond=None)
xycoefs,residsq,rank,s=np.linalg.lstsq(vec.T,R2,rcond=None)
rms=(residsq/len(xm))**0.5
R2fit=np.dot(vec.T,xycoefs)
ok=(R2fit-R2)**2<9*rms**2

print('clipping',(ok==False).sum(),'R2 values in fit over x,y')
# Iterate out outliers
xmfit=xm[ok]
ymfit=ym[ok]
one=np.ones(len(xmfit))
vec=np.vstack([np.ones(len(xmfit)),xmfit,ymfit,xmfit**2,xmfit*ymfit,ymfit**2,br[ok]])
xycoefs,residsq,rank,s=np.linalg.lstsq(vec.T,R2[ok],rcond=None)

vec=np.vstack([np.ones(len(xm)),xm,ym,xm**2,xm*ym,ym**2,br])
R2fit=np.dot(vec.T,xycoefs)-xycoefs[0]-xycoefs[-1]*br   # fit without the constant component & bprp, ie xy dependence only

plt.subplot(121)
plt.scatter(br,np.exp(frat),c=lf2,vmin=15.5,vmax=19.5,s=4,cmap='jet')
plt.plot(xx,np.exp(yy),'r')
plt.colorbar(label='Gaia G mag')
plt.xlabel('Gaia Bp-Rp')
plt.ylabel('VIS/Gaia flux')
plt.xlim(0.5,3)
plt.ylim(0.6,1.4)
plt.yscale('log')
plt.yticks(np.linspace(0.6,1.4,9))
plt.grid()
plt.title(momcatvis)
plt.text(0.55,1.35,'gain model '+led)

plt.subplot(222)
plt.hist(br,bins=12,range=[0.5,3])

plt.subplot(224)
plt.scatter(br,R2-R2fit,c=lf2,vmin=15.5,vmax=19.5,s=4,cmap='jet')
plt.xlabel('Gaia Bp-Rp')
plt.ylabel('R2')
plt.xlim(0.5,3)
plt.ylim(1.7,2.4)
plt.grid()

plt.savefig(momcatvis[:-8] + "_fluxrat_bprp.png")

spaces=str()
spaces=spaces.center(len(momcatvis)+len(gaiavis)+3)    # empty string of given length
spaces1='# cat1'+spaces[6:len(momcatvis)+1]+'  cat2'+spaces[7+len(gaiavis):]
spaces ='#'+spaces[1:]
print(spaces1+' keep   sigxy      dx        dy      ddx/dx     ddy/dx     ddx/dy     ddy/dy')
print(spaces +'   N     pix       pix       pix      -----------  pix/0.75 deg  ----------')
print(momcatvis,gaiavis,'%5i %7.3f %9.3f %9.3f %10.3f %10.3f %10.3f %10.3f' % (
    (keep.sum(),sigxy*36000)+
    tuple(coef[0,:]*36000)+
    tuple(coef[1:,:].flatten()*36000*0.75) ) )
print(momcatvis,'VIS/GaiaBpRp 0.75 1.75 2.75 %7.4f %7.4f %7.4f %9.2f' %
          (tuple(np.exp(fratref))+(texp,)))
