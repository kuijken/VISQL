import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from astropy.table import Table
import json
from matplotlib import use

if len(argv)==1:
    print("Usage: p3 plotzptFP.py file file file ... file label exptime")
    exit()

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'
use('Agg')

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

led='V'
gain=gainchoice[led]
refled='3'
refgain=gainchoice[refled]
gainfactor={e:gain[e]/refgain[e] for e in gain}

x=[]
y=[]
z=[]
G=[]
F=[]
br=[]

exptime=float(argv[-1])
label=argv[-2]

for frattab in argv[1:-2]:
    ftab=Table.read(frattab[:42]+'_frat.txt',format='ascii')
    ok=(ftab['TEXP']==exptime) & (ftab['G']<19.5)
    x+=list(ftab['XWORLD'][ok])
    y+=list(ftab['YWORLD'][ok])
    z+=list(ftab['FRAT'][ok]*[gainfactor[e] for e in ftab['EXT'][ok]])
    G+=list(ftab['G'][ok])
    F+=list(ftab['FLUX'][ok])
    br+=list(ftab['BpRp'][ok])

print('Found a total of',len(z),'Gaia/VIS stars for',label,'with exposure time',exptime)

# plt.plot(br,z,'k,')
# plt.show()

# frat=np.log(z)
# one=np.ones(len(br))
# vec=np.vstack([one,br,br**2])
# coefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([frat]).T,rcond=None)
# rms=(residsq/len(br))**0.5
# fit=coefs[0]+br*(coefs[1]+br*coefs[2])
# ok=(fit-frat)**2<9*rms**2
# print('clipping',(ok==False).sum(),'Gaia Bp-Rp stars')
# # Iterate out outliers
# brfit=br[ok]
# one=np.ones(len(brfit))
# vec=np.vstack([one,brfit,brfit**2])
# coefs,residsq,rank,s=np.linalg.lstsq(vec.T,np.vstack([frat[ok]]).T,rcond=None)
# fit=coefs[0]+br*(coefs[1]+br*coefs[2])
# z=np.exp(frat-fit)



def clip3n(aa):
    a=np.array(aa)
    ok=np.ones(len(a),dtype=bool)
    for iter in range(5):
        m=a[ok].mean()
        s=a[ok].std()
        ok=(a-m)**2<9*s**2
    return float(ok.sum())

def clip3m(aa):
    a=np.array(aa)
    ok=np.ones(len(a),dtype=bool)
    for iter in range(5):
        if ok.sum():
            m=a[ok].mean() 
            s=a[ok].std()
            ok=(a-m)**2<9*s**2
        else:
            m=np.nan
    return m

def clip3s100(aa):
    a=np.array(aa)
    ok=np.ones(len(a),dtype=bool)
    for iter in range(5):
        if ok.sum()>1:
            m=a[ok].mean()
            s=a[ok].std()
            ok=(a-m)**2<9*s**2
    if ok.sum()<2:
        s=np.nan
    return s*100

plt.figure(figsize=(13,8))
plt.subplot(231)
plt.hexbin(x,y,vmin=0,cmap='BuPu')
plt.colorbar(label='N')
plt.ylabel('Y_FPA')
plt.title('Data from '+label)

plt.subplot(232)
plt.hexbin(x,y,z,vmin=0.75,vmax=1.25,reduce_C_function=clip3m,cmap='jet')
plt.colorbar(label='mean VIS/Gaia pred')
plt.title('%i Gaia stars' % len(z))

plt.subplot(233)
plt.hexbin(x,y,z,vmin=0.98,vmax=1.02,reduce_C_function=clip3m,cmap='winter')
plt.colorbar(label='mean VIS/Gaia pred')
plt.title('Exposure time %6.2f' % exptime)

plt.subplot(234)
plt.hexbin(x,y,G,vmin=15,vmax=19.5,cmap='copper')
plt.colorbar(label='mean Gaia G')
plt.xlabel('X_FPA')
plt.ylabel('Y_FPA')

plt.subplot(235)
plt.hexbin(x,y,np.log10(F),vmin=4.4,vmax=5.4,cmap='copper')
plt.xlabel('X_FPA')
plt.colorbar(label='mean log10(VIS ADU)')

plt.subplot(236)
plt.hexbin(x,y,z,vmin=0,vmax=3,reduce_C_function=clip3s100,cmap='summer')
plt.colorbar(label='100 x sigma VIS/Gaia pred')
plt.title('Gain model '+led)
plt.xlabel('X_FPA')


plt.savefig('iceplot_FPA_'+label+('_%.2f'%exptime)+'.png')


