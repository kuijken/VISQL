from os import environ
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import use
from sys import argv
from glob import glob

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('Agg')

ccd=argv[1]
yyyymm=argv[2]
longshort=argv[3]
if longshort=='long':
    texp=560.52
else:
    texp=89.52

fratlist=glob('C_'+yyyymm+'*frat.txt')
fratlist.sort()
print('Found',len(fratlist),'*frat.txt files from',yyyymm,'.')

col=[]
rat=[]
tim=[]
for fratfile in fratlist:
    diso=fratfile[2:6]+'-'+fratfile[6:8]+'-'+fratfile[8:10]+'T'+fratfile[11:13]+':'+fratfile[13:15]+':'+fratfile[15:17]
    obsdate=datetime.datetime.fromisoformat(diso)
    obsdate=(obsdate - datetime.datetime(2024,1,1,0,0,0)).total_seconds() / (24*3600)
    z=Table.read(fratfile,format='ascii')
    onccd=z[(z['EXT']==ccd) & (z['G']<19.5) & (z['TEXP']==texp)]
    onccd.add_column([obsdate for x in range(len(onccd))],name='OBSDATE')
    #plt.scatter(onccd['BpRp'],onccd['VFLUX']/onccd['GFLUX'],c=onccd['OBSDATE'],s=3)
    col+=list(onccd['BpRp'])
    rat+=list(onccd['VFLUX']/onccd['GFLUX'])
    tim+=list(onccd['OBSDATE'])
#    print(obsdate)

plt.scatter(col,rat,c=tim,s=3)
plt.colorbar(label='days since 2024-01-01T0')
plt.title(ccd+'   '+yyyymm+'  '+longshort)
plt.xlabel('Bp-Rp')
plt.ylabel('VIS/G')
plt.savefig('frat-v-col.png')
plt.clf()

# fit fiducial colour eq to time 73.9 - 83 - period just after decontamination.

col=np.array(col)
rat=np.array(rat)
tim=np.array(tim)
f=(tim>73.9) & (tim<83)
vec=np.vstack([np.ones(f.sum()),col[f],col[f]**2])
tofit=np.log(rat[f])
ok=np.ones(f.sum(),dtype=bool)
for iter in range(4):
    coefs,residsq,rank,s=np.linalg.lstsq(vec.T[ok],tofit[ok],rcond=None)
    fit=coefs[0]+col[f]*(coefs[1]+coefs[2]*col[f])
    rmsq=residsq[0]/ok.sum()
    ok=(fit-tofit)**2<9*rmsq
    # print('clipping',(ok==False).sum(),'Gaia Bp-Rp stars',rmsq**0.5,ok.sum())
plt.plot(col[f],tofit,'k.')
plt.plot(col[f],fit,'r.')
plt.savefig('colcolfit.png')
plt.clf()


c0,c1,c2= -0.445329 ,  0.402658 , -0.063669
c0,c1,c2=coefs
print(c0,c1,c2)
fiducial=np.exp(c0+col*(c1+col*c2))
plt.figure(figsize=(10,6))
plt.scatter(tim,2.5*np.log10(rat/fiducial),c=col,s=3,vmin=0.5,vmax=3,cmap='jet')
plt.colorbar(label='Bp-Rp')
plt.grid()
plt.xlabel('Days since 2024-01-01T00')
plt.ylabel('VIS zpt vs March 14-22 2024 [mag]')
plt.title(ccd+'   '+yyyymm+'  '+longshort)
plt.ylim(-0.7,0.1)
plt.savefig('iceplot_%s_%s_%.2f.png' % (ccd,yyyymm,texp))
plt.clf()
