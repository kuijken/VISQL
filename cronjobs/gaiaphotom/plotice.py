from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import use
import datetime
import numpy as np

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=9
plt.rcParams['ytick.labelsize']=9
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'
use('Agg')

t=Table.read('long.txt',format='ascii')

# generate dates to plot on x axis
daytab=np.array(
    [datetime.datetime.fromisoformat(
        x[2:6]+'-'+x[6:8]+'-'+x[8:10]+'T'+x[11:13]+':'+x[13:15]+':'+x[15:17])
    for x in t['col1']
    ]
    )

plt.figure(figsize=(8,6))
plt.title('Long exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(-0.4,0.2)
plt.xlim(datetime.date(2023,7,15),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_long.png', bbox_inches="tight")
plt.clf()

plt.figure(figsize=(6,6))
plt.title('Long exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(-0.4,0.2)
plt.xlim(datetime.date(2024,3,1),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_long2.png', bbox_inches="tight")
plt.clf()

plt.figure(figsize=(6,6))
plt.title('Long exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(0,0.2)
plt.xlim(datetime.date(2024,6,1),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_long3.png', bbox_inches="tight")
plt.clf()


t=Table.read('short.txt',format='ascii')
daytab=np.array(
    [datetime.datetime.fromisoformat(
        x[2:6]+'-'+x[6:8]+'-'+x[8:10]+'T'+x[11:13]+':'+x[13:15]+':'+x[15:17])
    for x in t['col1']
    ]
    )
plt.figure(figsize=(8,6))
plt.title('Short exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(-0.4,0.2)
plt.xlim(datetime.date(2023,7,15),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_short.png', bbox_inches="tight")
plt.clf()

plt.figure(figsize=(6,6))
plt.title('Short exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(-0.4,0.2)
plt.xlim(datetime.date(2024,3,1),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_short2.png', bbox_inches="tight")
plt.clf()

plt.figure(figsize=(6,6))
plt.title('Short exposures')
plt.scatter(daytab,np.log(t['col8']*0.95),c='r',label='Bp-Rp=2.75',s=0.1)
plt.scatter(daytab,np.log(t['col7']*1.),c='g',label='Bp-Rp=1.75',s=0.1)
plt.scatter(daytab,np.log(t['col6']*1.2),c='b',label='Bp-Rp=0.75',s=0.1)
plt.ylim(0,0.2)
plt.xlim(datetime.date(2024,6,1),datetime.date.today()+datetime.timedelta(days=7))
plt.xticks(rotation=75)
plt.grid()
plt.legend()
plt.tight_layout()
plt.ylabel('ln (VIS/Gaia)')
plt.savefig('iceplot_short3.png', bbox_inches="tight")
plt.clf()
