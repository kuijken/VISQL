from astropy.table import Table
import astropy.io.ascii as asc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
import datetime
from glob import glob
from sys import argv

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'
use('agg')

yyyymm=argv[1]

f=Table.read('psfmoments_'+yyyymm+'.txt',format='ascii.commented_header')

exptimes=asc.read('../FITS/fitslist.txt')
imname=f['catalogue']
texp=[exptimes['EXPTIME'][exptimes['FITSNAME']==im[:-15]+".fits"][0] for im in imname]

missionday=np.array([
    (datetime.datetime(
        int(n[2:6]),int(n[6:8]),int(n[8:10]),
        int(n[11:13]),int(n[13:15]),int(n[15:17])) -
        datetime.datetime(2023,7,1,0,0,0)).total_seconds()/86400.
    for n in imname])
day0=int(missionday[0])
dayf=int(missionday[-1])+1
dayf=day0+31
year0=int(imname[0][2:6])
month0=int(imname[0][6:8])
if month0 < 12:
    yearf=year0+0
    monthf=month0+1
else:
    yearf=year0+1
    monthf=1
day0=(datetime.date(year0,month0,1)-datetime.date(2023,7,1)      ).total_seconds()/86400.
dayf=(datetime.date(yearf,monthf,1)-datetime.date(2023,7,1)      ).total_seconds()/86400.
def mm(v,v1,v2):
    return np.minimum(np.maximum(v,v1),v2)

plt.figure(figsize=(14,12))
plt.subplot(511)
plt.scatter(missionday,mm(f['R2'],1.72,2.28),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(1.7,2.3)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('R2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(512)
plt.scatter(missionday,mm(f['ELL 1'],-0.065,0.065),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.07,0.07)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('e1')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(513)
plt.scatter(missionday,mm(f['ELL 2'],-0.065,0.065),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.07,0.07)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('e2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(514)
plt.scatter(missionday,mm(f['TREF 1'],-0.29,-0.16),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.3,-0.15)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('t1')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(515)
plt.scatter(missionday,mm(f['TREF 2'],-0.04,0.09),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.05,0.1)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('t2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.xlabel('Mission day ('+imname[0][2:17]+' to '+imname[-1][2:17]+')')
plt.tight_layout()
plt.savefig('IQ_timeline_'+yyyymm+'.png')
plt.clf()

# now repeat the plot, but plot last 31 days instead of current month as above

# first read in the previous month

y,m=int(yyyymm[:4]),int(yyyymm[4:6])
if (m==1):
    y-=1
    m=12
else:
    m-=1
yyyymm2='%04i%02d' % (y,m)
f2=Table.read('psfmoments_'+yyyymm2+'.txt',format='ascii.commented_header')

imname2=f2['catalogue']
texp2=[exptimes['EXPTIME'][exptimes['FITSNAME']==im[:-15]+".fits"][0] for im in imname2]

missionday2=np.array([
    (datetime.datetime(
        int(n[2:6]),int(n[6:8]),int(n[8:10]),
        int(n[11:13]),int(n[13:15]),int(n[15:17])) -
        datetime.datetime(2023,7,1,0,0,0)).total_seconds()/86400.
    for n in imname2])

dayf=int(missionday[-1])+1
day0=dayf-31
date0=imname[-1][2:10]+'_000000'
y,m,d=int(date0[:4]),int(date0[4:6]),int(date0[6:8])
date0=datetime.date(y,m,d)-datetime.timedelta(days=30)
date0='%04d%02d%02d_000000' % (date0.year,date0.month,date0.day)

plt.figure(figsize=(14,12))
plt.subplot(511)
plt.scatter(missionday,mm(f['R2'],1.72,2.28),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.scatter(missionday2,mm(f2['R2'],1.72,2.28),c=texp2,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(1.7,2.3)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('R2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(512)
plt.scatter(missionday,mm(f['ELL 1'],-0.065,0.065),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.scatter(missionday2,mm(f2['ELL 1'],-0.065,0.065),c=texp2,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.07,0.07)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('e1')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(513)
plt.scatter(missionday,mm(f['ELL 2'],-0.065,0.065),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.scatter(missionday2,mm(f2['ELL 2'],-0.065,0.065),c=texp2,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.07,0.07)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('e2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(514)
plt.scatter(missionday,mm(f['TREF 1'],-0.29,-0.16),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.scatter(missionday2,mm(f2['TREF 1'],-0.29,-0.16),c=texp2,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.3,-0.15)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('t1')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.subplot(515)
plt.scatter(missionday,mm(f['TREF 2'],-0.04,0.09),c=texp,vmin=0,vmax=700,s=1,cmap='turbo')
plt.scatter(missionday2,mm(f2['TREF 2'],-0.04,0.09),c=texp2,vmin=0,vmax=700,s=1,cmap='turbo')
plt.ylim(-0.05,0.1)
plt.xlim(day0,dayf)
plt.xticks(np.arange(day0,dayf+1,1))
plt.ylabel('t2')
plt.grid()
plt.colorbar(label='Texp',aspect=10)
plt.xlabel('Mission day ('+date0+' to '+imname[-1][2:17]+')')
plt.tight_layout()
plt.savefig('IQ_timeline_lastmonth.png')
plt.clf()

