import astropy.io.ascii as asc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
import datetime
from glob import glob
from sys import argv
from astropy.table import Table
from datetime import date

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=9
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'
use('agg')

if len(argv)<2:
    print("Usage: p3 plotfratdiagperday.py day1 day2 long/short")
    print("   where days are in format yyyymmdd")
    print(" This will read in vis-over-gaia-yyyymmdd-ccd_long/short.txt files and plot")
    print("   the flux loss in the diagonal quadrants between day 1 and day2")
    exit()
    
day1,day2=argv[1:3]
y1=int(day1[:4])
m1=int(day1[4:6])
d1=int(day1[6:8])
y2=int(day2[:4])
m2=int(day2[4:6])
d2=int(day2[6:8])
longshort=argv[3]
# first set reference as avg 20240316-21
dplot=[]
rat1=[]
rat2=[]
rat3=[]
rat4=[]
rat5=[]
rat6=[]
xrat1=[]
xrat2=[]
xrat3=[]
xrat4=[]
xrat5=[]
xrat6=[]
for d in range(6):
  try:
    day=date(2024,3,16)+datetime.timedelta(days=d)
    daystring='%04i%02i%02i' % (day.year,day.month,day.day)
    t=Table.read('vis-over-gaia-'+daystring+'-ccd_'+longshort+'.txt',format='ascii')
#    quads=[str(n)+'-'+str(n)+'.GVIS/Gaia' for n in range(1,7)]
#    rat=[t['VoverG2'][t['RATIONAME']==q] for q in quads]
    ccds=['CCD'+str(n)+str(n)+'VIS/Gaia' for n in range(1,7)]
    rat=[t['VoverG2'][t['RATIONAME']==c] for c in ccds]
    dplot+=[day]
    rat1+=[rat[0]]
    rat2+=[rat[1]]
    rat3+=[rat[2]]
    rat4+=[rat[3]]
    rat5+=[rat[4]]
    rat6+=[rat[5]]
    ccds=['CCD'+str(n)+str(7-n)+'VIS/Gaia' for n in range(1,7)] # other diagonal
    rat=[t['VoverG2'][t['RATIONAME']==c] for c in ccds]
    xrat1+=[rat[0]]
    xrat2+=[rat[1]]
    xrat3+=[rat[2]]
    xrat4+=[rat[3]]
    xrat5+=[rat[4]]
    xrat6+=[rat[5]]
  except:
    pass

ref1=np.array(rat1).mean()
ref2=np.array(rat2).mean()
ref3=np.array(rat3).mean()
ref4=np.array(rat4).mean()
ref5=np.array(rat5).mean()
ref6=np.array(rat6).mean()
xref1=np.array(xrat1).mean()
xref2=np.array(xrat2).mean()
xref3=np.array(xrat3).mean()
xref4=np.array(xrat4).mean()
xref5=np.array(xrat5).mean()
xref6=np.array(xrat6).mean()

# now get the data to be plotted
dplot=[]
rat1=[]
rat2=[]
rat3=[]
rat4=[]
rat5=[]
rat6=[]
xrat1=[]
xrat2=[]
xrat3=[]
xrat4=[]
xrat5=[]
xrat6=[]
for d in range((date(y2,m2,d2)-date(y1,m1,d1)).days+1):
  try:
    day=date(y1,m1,d1)+datetime.timedelta(days=d)
    daystring='%04i%02i%02i' % (day.year,day.month,day.day)
    t=Table.read('vis-over-gaia-'+daystring+'-ccd_'+longshort+'.txt',format='ascii')
#    quads=[str(n)+'-'+str(n)+'.GVIS/Gaia' for n in range(1,7)]
#    rat=[t['VoverG2'][t['RATIONAME']==q] for q in quads]
    ccds=['CCD'+str(n)+str(n)+'VIS/Gaia' for n in range(1,7)]
    rat=[t['VoverG2'][t['RATIONAME']==c] for c in ccds]
    dplot+=[day]
    rat1+=[rat[0]/ref1]
    rat2+=[rat[1]/ref2]
    rat3+=[rat[2]/ref3]
    rat4+=[rat[3]/ref4]
    rat5+=[rat[4]/ref5]
    rat6+=[rat[5]/ref6]
    ccds=['CCD'+str(n)+str(7-n)+'VIS/Gaia' for n in range(1,7)]
    xrat=[t['VoverG2'][t['RATIONAME']==c] for c in ccds]
    xrat1+=[xrat[0]/xref1]
    xrat2+=[xrat[1]/xref2]
    xrat3+=[xrat[2]/xref3]
    xrat4+=[xrat[3]/xref4]
    xrat5+=[xrat[4]/xref5]
    xrat6+=[xrat[5]/xref6]
  except:
    pass

plt.subplot(2,1,1)
plt.plot(dplot,rat1,label='1-1')
plt.plot(dplot,rat2,label='2-2')
plt.plot(dplot,rat3,label='3-3')
plt.plot(dplot,rat4,label='4-4')
plt.plot(dplot,rat5,label='5-5')
plt.plot(dplot,rat6,label='6-6')
#plt.legend()
plt.xticks(labels=None,rotation=12)
plt.grid()
#plt.xlabel('date')
plt.ylabel('VIS/Gaia')
#plt.title('diagonal G quadrants - '+longshort)
plt.title('diagonal CCDs - '+longshort)
plt.subplot(2,1,2)
plt.plot(dplot,rat1,label='1-1')
plt.plot(dplot,rat2,label='2-2')
plt.plot(dplot,rat3,label='3-3')
plt.plot(dplot,rat4,label='4-4')
plt.plot(dplot,rat5,label='5-5')
plt.plot(dplot,rat6,label='6-6')
plt.ylim(0.9,1.05)
plt.legend()
plt.xticks(rotation=12)
plt.grid()
plt.xlabel('date')
plt.ylabel('VIS/Gaia')
#plt.title('diagonal G quadrants - '+longshort)
#plt.title('diagonal CCDs - '+longshort)
plt.tight_layout()
plt.savefig('fratdiagtimeline.png')
plt.clf()

plt.subplot(2,1,1)
plt.plot(dplot,xrat1,label='1-6')
plt.plot(dplot,xrat2,label='2-5')
plt.plot(dplot,xrat3,label='3-4')
plt.plot(dplot,xrat4,label='4-3')
plt.plot(dplot,xrat5,label='5-2')
plt.plot(dplot,xrat6,label='6-1')
#plt.legend()
plt.xticks(labels=None,rotation=12)
plt.grid()
#plt.xlabel('date')
plt.ylabel('VIS/Gaia')
#plt.title('cross diagonal G quadrants - '+longshort)
plt.title('cross diagonal CCDs - '+longshort)
plt.subplot(2,1,2)
plt.plot(dplot,xrat1,label='1-6')
plt.plot(dplot,xrat2,label='2-5')
plt.plot(dplot,xrat3,label='3-4')
plt.plot(dplot,xrat4,label='4-3')
plt.plot(dplot,xrat5,label='5-2')
plt.plot(dplot,xrat6,label='6-1')
plt.ylim(0.9,1.05)
plt.legend()
plt.xticks(rotation=12)
plt.grid()
plt.xlabel('date')
plt.ylabel('VIS/Gaia')
#plt.title('diagonal G quadrants - '+longshort)
#plt.title('diagonal CCDs - '+longshort)
plt.tight_layout()
plt.savefig('fratxdiagtimeline.png')
plt.clf()
