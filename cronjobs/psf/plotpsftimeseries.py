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


