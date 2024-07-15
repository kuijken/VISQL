import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
from sys import argv
import datetime
import os, requests
import netCDF4 as nc
from astropy.time import Time
import time
import astropy.io.fits as pf
import os
from os import environ


# get goes data for a single yyy mm dd, save to txt file


t0euclid=datetime.datetime(2023,7,1,0,0,0)
t0goes  =datetime.datetime(2000,1,1,12,0,0)
timezeropoint=(t0euclid-t0goes).total_seconds()


def get_goes_data(year,month,day,verbose=False):
        
    dat = datetime.datetime(year,month,day).strftime('%Y%m%d')
    filename = 'sci_xrsf-l2-avg1m_g16_d'+dat+'_v2-2-0.nc'
    asciifile =environ['viskom']+'/XRAY/goes'+dat+'.txt'
    #only download if you haven't already
    t,flux=[],[]
    if not os.path.exists(asciifile): 
            with open(filename, "wb") as f:
                url_path = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-avg1m_science/"+dat[0:4]+'/'+dat[4:6]+'/'
                r = requests.get(url_path+filename)
                if r.reason=='Not Found': #if the data is missing, print the name and delete the empty file
                    print('Could not download file: '+filename)
                    open(asciifile,'w').close()
                else:
                    f.write(r.content)
                    f.close()
                    ff=nc.Dataset(filename)
                    t=ff.variables['time'][:].data           #time in seconds since 2000-01-01 12:00:00
                    flux = ff.variables['xrsa_flux'][:].data #hard X-ray flux
                    flag = ff.variables['xrsa_flag'][:].data # flag
                    t=t[flag==0]-timezeropoint
                    flux=flux[flag==0]
                    np.savetxt(asciifile,np.array([t,flux]).T,fmt='%10.f %10.5g',
                                header='time since 2023-07-01T00:00:00 (s);   GOES XRS-A FLUX (W/m2)')
                    if verbose:
                        print('Read GOES data from NOAA and saved in file %s (%i values)' % (asciifile,len(t)))
                os.remove(filename)
    else:
        try:
            t,flux=np.loadtxt(asciifile,unpack=True)
            if verbose:
                print('Read GOES data from existing file %s (%i values)' % (asciifile,len(t)))
        except:
            t,flux=[],[]
            if verbose:
                print('Read GOES data from existing file %s (%i values)' % (asciifile,len(t)))
            
    return t,flux

def get_goes_data_raw(year,month,day,verbose=False):
        
    dat = datetime.datetime(year,month,day).strftime('%Y%m%d')
    filename = 'dn_xrsf-l2-avg1m_g16_d'+dat+'_v2-2-0.nc'
    asciifile =environ['viskom']+'/XRAY/goes'+dat+'.txt'
    #only download if you haven't already
    t,flux=[],[]
    if not os.path.exists(asciifile): 
            with open(filename, "wb") as f:
                url_path = "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-avg1m/"+dat[0:4]+'/'+dat[4:6]+'/'
                r = requests.get(url_path+filename)
                if r.reason=='Not Found': #if the data is missing, print the name and delete the empty file
                    print('Could not download file: '+filename)
                    open(asciifile,'w').close()
                else:
                    f.write(r.content)
                    f.close()
                    ff=nc.Dataset(filename)
                    t=ff.variables['time'][:].data           #time in seconds since 2000-01-01 12:00:00
                    flux = ff.variables['xrsa_flux'][:].data #hard X-ray flux
                    flag = ff.variables['xrsa_flag'][:].data # flag
                    flag*=0
                    t=t[flag==0]-timezeropoint
                    flux=flux[flag==0]
                    np.savetxt(asciifile,np.array([t,flux]).T,fmt='%10.f %10.5g',
                                header='time since 2023-07-01T00:00:00 (s);   GOES XRS-A FLUX (W/m2)')
                    if verbose:
                        print('Read GOES data from NOAA and saved in file %s (%i values)' % (asciifile,len(t)))
                os.remove(filename)
    else:
        try:
            t,flux=np.loadtxt(asciifile,unpack=True)
            if verbose:
                print('Read GOES data from existing file %s (%i values)' % (asciifile,len(t)))
        except:
            t,flux=[],[]
            if verbose:
                print('Read GOES data from existing file %s (%i values)' % (asciifile,len(t)))
            
    return t,flux

t,f=get_goes_data_raw(int(argv[1]),int(argv[2]),int(argv[3]),verbose=True)

