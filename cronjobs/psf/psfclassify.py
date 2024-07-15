from sys import argv
import matplotlib.pyplot as plt
from matplotlib import use
import numpy as np
from os import environ

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')
#use('MacOSX')

'''
Read in a table listing median PSF moments for a set of observations.
(produced with medianmoments).
Then look for outliers in the various aberrations
'''

nsig=2.5

mediantablefile=argv[1]
r2,c1,c2,e1,e2,t1,t2=np.loadtxt(mediantablefile,unpack=True,usecols=[3,5,7,9,11,13,15])
im=np.loadtxt(mediantablefile,unpack=True,usecols=[1],dtype=str)

coma=np.array([c1,c2]).T
ell= np.array([e1,e2]).T
tref=np.array([t1,t2]).T
# stats for r2
r2ok = r2<2.5
for iter in range(10):
    r2rms=np.std(r2[r2ok])
    r2mean=np.mean(r2[r2ok])
    r2ok=np.abs(r2-r2mean)<nsig*r2rms
    print('R2',iter,r2mean,r2rms,r2ok.sum())
# stats for coma (do this in 2d)
comaok= c1**2+c2**2 < 0.03**2
for iter in range(10):
    comamean=coma[comaok].mean(axis=0)
    comavar=coma[comaok].var(axis=0).sum()
    comaok=((coma-comamean)**2).sum(axis=1)<comavar*nsig**2
    print('COMA',iter,comamean,comavar**0.5,comaok.sum())
    
ellok= e1**2+e2**2 < 0.07**2
for iter in range(10):
    ellmean=ell[ellok].mean(axis=0)
    ellvar=ell[ellok].var(axis=0).sum()
    ellok=((ell-ellmean)**2).sum(axis=1)<ellvar*nsig**2
    print('ELL',iter,ellmean,ellvar**0.5,ellok.sum())

trefok= t1**2+t2**2 < 0.3**2
for iter in range(10):
    trefmean=tref[trefok].mean(axis=0)
    trefvar=tref[trefok].var(axis=0).sum()
    trefok=((tref-trefmean)**2).sum(axis=1)<trefvar*nsig**2
    print('TREF',iter,trefmean,trefvar**0.5,trefok.sum())
    
print((r2ok & comaok & ellok & trefok).sum())

oksumm=(r2<1.95) & (r2>1.6) & ((e1+0.02)**2+e2**2<0.012**2)
score= (r2-1.9)**2/0.1**2 + (e1+0.02)**2/0.007**2 + e2**2/0.007**2
score= (r2-1.875)**2/0.025**2 + (e1+0.02)**2/0.01**2 + e2**2/0.01**2
#score= np.minimum(r2-1.8,np.minimum(0,1.9-r2))**2/0.025**2  + (e1+0.017)**2/0.01**2 + e2**2/0.01**2

ok= r2<3
if ok.sum()>10:
    for iter in range(10):
        r2m=r2[ok].mean()
        e1m=e1[ok].mean()
        e2m=e2[ok].mean()
        vare2=e2[ok].var()
        cov=np.cov(r2[ok],e1[ok])
        covinv=np.linalg.inv(cov)
        print(r2m,e1m,e2m,ok.sum(),cov[0,0],cov[1,0],cov[1,1],vare2)
        ok=(e2-e2m)**2/vare2 + (r2-r2m)**2*covinv[0,0]+2*(r2-r2m)*(e1-e1m)*covinv[1,0]+(e1-e1m)**2*covinv[1,1]<3*2.5**2
    print(r2m,e1m,e2m,ok.sum(),cov[0,0],cov[1,0],cov[1,1],vare2)


# calibrated on 20230920-20230923_230344:

r2m=1.9233812949640285
e1m=-0.023129496402877696
e2m=0.002733812949640288
cov00=0.0005422666041080192
cov01=-0.00011841403399019925
cov11=4.4881659889479756e-05
vare2=5.332022152062523e-06
print(r2m,e1m,e2m,ok.sum(),cov00,cov01,cov11,vare2,'adopted 20230923')

covinv=np.linalg.inv([[cov00,cov01],[cov01,cov11]])

score=(e2-e2m)**2/vare2 + (r2-r2m)**2*covinv[0,0]+2*(r2-r2m)*(e1-e1m)*covinv[1,0]+(e1-e1m)**2*covinv[1,1]

oksumm=score<27

plt.subplot(221)
plt.plot(r2,e1,'b,')
plt.plot(r2,e2,'r,')
plt.plot(r2[oksumm],e1[oksumm],'b.',label='e1')
plt.plot(r2[oksumm],e2[oksumm],'r.',label='e2')
plt.xlim(1.8,2.1)
plt.ylim(-0.05,0.025)
plt.xlabel('R2')
plt.ylabel('ellip1,2')
plt.legend()

plt.subplot(222)
plt.plot(e1,e2,'k.')
plt.plot(e1[oksumm],e2[oksumm],'r.')
plt.xlabel('ellip1')
plt.ylabel('ellip2')
plt.xlim(-0.05,0.0)
plt.ylim(-0.025,0.025)

plt.subplot(223)
plt.plot(c1,c2,'k.')
plt.plot(c1[oksumm],c2[oksumm],'r.')
plt.xlabel('coma1')
plt.ylabel('coma2')
plt.xlim(0,0.05)
plt.ylim(-0.05,0)

plt.subplot(224)
plt.plot(t1,t2,'k.')
plt.plot(t1[oksumm],t2[oksumm],'r.')
plt.xlabel('tref1')
plt.ylabel('tref2')
plt.xlim(-0.3,-0.2)
plt.ylim(-0.05,0.05)

plt.tight_layout()
plt.savefig(mediantablefile[:-4]+".png")
plt.close()

outfilename=mediantablefile[:-4]+"_score.txt"
outfile=open(outfilename,'w')
print('| image | R2 | e1 | e2 | score |',file=outfile)
for i in range(len(im)):
    print('|',im[i][:-14],' | %8.3f | %8.4f | %8.4f | %8.2f |' % (r2[i],e1[i],e2[i],score[i]),file=outfile)
outfile.close()

csvfilename=mediantablefile[:-4]+"_score.csv"
csvfile=open(csvfilename,'w')
print('image , R2 , e1 , e2 , score ',file=csvfile)
for i in range(len(im)):
    print(im[i][:-14],' , %8.3f , %8.4f , %8.4f , %8.2f' % (r2[i],e1[i],e2[i],score[i]),file=csvfile)
csvfile.close()

le1name,iwsname=np.loadtxt(environ['viskom']+'/FITS/le1-vs-iws.txt',unpack=True,dtype=str)

le1=[(list(le1name[iwsname==image[:-15]+".fits"])+[""])[0] for image in im]

outfilename=mediantablefile[:-4]+"_le1_score.txt"
outfile=open(outfilename,'w')
print('| image | R2 | e1 | e2 | score | EASname |',file=outfile)
for i in range(len(im)):
    print('|',im[i][:-14],' | %8.3f | %8.4f | %8.4f | %8.2f |' % (r2[i],e1[i],e2[i],score[i]),le1[i],' | ',file=outfile)
outfile.close()

csvfilename=mediantablefile[:-4]+"_le1_score.csv"
csvfile=open(csvfilename,'w')
print('image ; R2 ; e1 ; e2 ; score ; EASname ',file=csvfile)
for i in range(len(im)):
    print(im[i][:-14],' ; %8.3f ; %8.4f ; %8.4f ; %8.2f' % (r2[i],e1[i],e2[i],score[i]),"; ",le1[i],file=csvfile)
csvfile.close()


