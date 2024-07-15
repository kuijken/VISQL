import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
import os 
from sys import argv
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.labelweight']='normal'

use('agg')

t0=Time('2023-07-01T00:00:00')
readtime=72 # readout time
viskom=os.environ['viskom']
angxtab=Table.read(viskom+'/XRAY/All_Angles_and_GOES.txt',format='ascii')
tx,xgoes=np.loadtxt(viskom+'/XRAY/AllGOES.txt',unpack=True)

for i in range(len(angxtab)):
    pngfile=angxtab['FITSNAME'][i][:-5]+'_goes.png'
    if os.path.exists(pngfile):
        #print(pngfile,'already exists.')
        pass
    else:
        dobs=Time(angxtab['DATE-OBS'][i])
        timewindow=(t0+tx*u.s>dobs-12*u.h) & (t0+tx*u.s<dobs+12*u.h)
        plt.plot(tx[timewindow]/3600+(t0-dobs).to_value('h'),xgoes[timewindow])
        plt.plot([0,0],[0,1],'r')
        plt.xlim(-12,12)
        plt.ylim(1e-8,1e-4)
        plt.xlabel('Hours from observation')
        plt.ylabel('GOES flux W/m2')
        plt.yscale('log')
        plt.title(angxtab['FITSNAME'][i][:-5])
        plt.text(10,3e-5,'%5.1f sec + %4.0f' %(angxtab['EXPTIME'][i],readtime),ha='right')
        plt.text(-11,5e-5,'SAA  =%7.2f' % angxtab['SAA'][i])
        plt.text(-11,3e-5,'ALPHA=%7.2f' % angxtab['ALPHA'][i])
        plt.savefig(pngfile)
        plt.close()
        print('Written',pngfile)
