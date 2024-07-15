from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import numpy as np
import datetime as dt
from sys import argv
from os import environ

time0mjd=2460126.5  # time zero point in Herve's table (MJD) - Euclid time (days from 2023-07-01T00:00:00)
#dt0=dt.datetime(2023,7,1,0,0,0)

# Need to modify this to look at header of fits files and the SAA-ALFA-recon table.
# TBD: Include le1 header information as first option!

readout=72

yyyymm=argv[1]

GOESname  =environ['viskom']+'/XRAY/AllGOES.txt'
anglename =environ['viskom']+'/XRAY/SAA-ALFA-recon.txt'
le1angname=environ['viskom']+'/XRAY/SAA-ALFA-LE1.txt'
fitsname  =environ['viskom']+'/FITS/fitslist.txt'
srawname  =environ['viskom']+'/XRAY/sraw_2023_08_22_with_times.fits'
pvname    =environ['viskom']+'/XRAY/pv_list_from_20230812T000000_to_20230818T094336.fits'

fitstab =Table.read(fitsname,format='ascii')
angletab=Table.read(anglename,format='ascii')
le1tab  =Table.read(le1angname,format='ascii')
goestab =Table.read(GOESname,format='ascii')
srawtab =Table.read(srawname)
pvtab   =Table.read(pvname)

goestime=goestab['t_euclid_sec']*u.s+Time('2023-07-01')
print('Latest GOES time:',goestime[-1])
print('Latest angle    :',angletab['filename'][-1])
print('Latest LE1 angle:',le1tab['filename'][-1])

srawfile=srawtab['SRAW File']
srawiexp=srawtab['iExp']
srawfits=np.array([srawfile[i]+(".bin_%02d_01.fits" %srawiexp[i]) for i in range(len(srawfile))])

pvfile=pvtab['FileName']
pvfits1=np.array([ p[20:30]+"_"+p[31:37] for p in pvfile])
pvfits2=np.array([ p[45:51]+".fits"  for p in pvfile])

# find fits files in commandline-specified month
inyyyymm=[f[2:8]==yyyymm for f in fitstab['FITSNAME']]
fitstab=fitstab[inyyyymm]
print('Processing',len(fitstab),'files from',yyyymm)

texp=fitstab['EXPTIME']
filename=fitstab['FITSNAME']
Nfits=len(filename)
dateobs=fitstab['DATE-OBS']
goespersec=np.zeros(Nfits)
goestotal =np.zeros(Nfits)
saa       =np.zeros(Nfits)
alpha     =np.zeros(Nfits)
ra        =np.zeros(Nfits)
dec       =np.zeros(Nfits)
pa        =np.zeros(Nfits)
source    =np.zeros(Nfits,dtype='<U8')
for i in range(Nfits):
    fitsfile=filename[i]
    dstart=Time(dateobs[i])
    if dstart < Time('2000-01-01'):
        timefromfilename=fitsfile[2:6]+"-"+fitsfile[6:8]+"-"+fitsfile[8:10]+ \
            "T"+fitsfile[11:13]+":"+fitsfile[13:15]+":"+fitsfile[15:17]
        dstart=Time(timefromfilename)
        dateobs[i]=dstart.isot
    dend=dstart+(readout+texp[i])*u.s
    # find the GOES data during the exposure.
    ingoes=(goestime>=dstart) & (goestime<=dend)
    if ingoes.sum():
        goespersec[i]=goestab['GOES_X'][ingoes].mean()
    else:
        goespersec[i]=goestab['GOES_X'][np.argmin(np.abs(dstart-goestime))]
    goestotal[i]=goespersec[i]*(readout/2+texp[i])
    # now look for angles. These are labelled by .png thumnail filenames.
    # first look for LE1 header information
    inangle=le1tab['filename']==fitsfile[:-5]+".png"
    inangle2=angletab['filename']==fitsfile[:-5]+".png"
    print(fitsfile)
    if inangle.sum():
        iang=np.argmax(inangle)
        source[i]='LE1'
        saa[i]=le1tab['SAA'][iang]
        alpha[i]=le1tab['ALPHA'][iang]
        ra[i]=le1tab['RA'][iang]        
        dec[i]=le1tab['DEC'][iang]       
        pa[i]=le1tab['PA'][iang]
        print(filename[i],'has LE1 angle info: index',iang)
    else:
        print('No match for',fitsfile)
        source[i]='NONE'
        saa[i]=999
        alpha[i]=999
        ra[i]=999
        dec[i]=999
        pa[i]=999

print('Done trying to obtain angles for all files')

# fitstab['SAA']  =saa
# fitstab['ALPHA']=alpha     
# fitstab['RA']   =ra        
# fitstab['DEC']  =dec       
# fitstab['PA']   =pa        
# fitstab['Xgoes_per_sec']   =goespersec
# fitstab['Xgoes_integrated']=goestotal
# fitstab['ANGSRC']=source

#out=open(environ['viskom']+'/XRAY/All_Angles_and_GOES.txt','w')

#print("FITSNAME DATE-OBS EXPTIME IMG_T1 IMG_T2 SAA ALPHA RA DEC PA Xgoes_per_sec Xgoes_integrated ANGSRC",file=out)

#for i in range(len(filename)):
#    print('%s %s %8.2f %s %s %12.7f %11.7f %12.7f %12.7f %12.7f %12.5g %12.5g %s' %
#          (fitstab['FITSNAME'][i],fitstab['DATE-OBS'][i],fitstab['EXPTIME'][i],
#           fitstab['IMG_T1'][i],fitstab['IMG_T2'][i],
#           saa[i],alpha[i],ra[i],dec[i],pa[i],goespersec[i],goestotal[i],source[i]),file=out)
#out.close()
#fitstab.write(environ['viskom']+'/XRAY/All_Angles_and_GOES.txt',format='ascii',overwrite=True)

out=Table([fitstab['FITSNAME'],fitstab['DATE-OBS'],fitstab['EXPTIME'],fitstab['IMG_T1'],fitstab['IMG_T2']])
out.add_column(saa,name ='SAA')
out.add_column(alpha,name='ALPHA')
out.add_column(ra,name='RA')
out.add_column(dec,name='DEC')
out.add_column(pa,name='PA')
out.add_column(goespersec,name='Xgoes_per_sec')
out.add_column(goestotal,name='Xgoes_integrated')
out.add_column(source,name='ANGSRC')

out.write(environ['viskom']+'/XRAY/All_Angles_and_GOES_'+yyyymm+'.txt',overwrite=True,format='ascii')

