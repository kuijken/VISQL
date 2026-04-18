import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import use
from sys import argv
import datetime

# select X rays as having flux between 400 and 600 ADU
# plot counts on VIS FP

plt.rcParams['text.usetex']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['axes.titlesize']=10
plt.rcParams['axes.labelweight']='normal'

use('agg')

t=Table.read(argv[1],format='ascii.no_header',names=['time','filename'])
title=argv[2]
times=np.array([datetime.datetime.fromisoformat(x) for x in t['time']])

# locate specific file types and plot their creation time as a stepped line

keep = np.strings.find(t['filename'],'.fits') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='FITS')
print('FITS',keep.sum())

keep = np.strings.find(t['filename'],'01.png') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='THUMB')
print('THUMB',keep.sum())

keep = np.strings.find(t['filename'],'.cat.gz') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='CAT')
print('CAT',keep.sum())

keep = np.strings.find(t['filename'],'mom.cat') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='MOMCAT')
print('MOMCAT',keep.sum())

keep = np.strings.find(t['filename'],'wcs') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='WCS')
print('WCS',keep.sum())

keep = np.strings.find(t['filename'],'iq.png') >-1
i=np.arange(keep.sum())+1
plt.plot(times[keep],i,ds='steps-pre',label='IQ')
print('IQ',keep.sum())

plt.xlabel('File creation time')
plt.ylabel('Cumul')
plt.title(title)
plt.legend()
plt.xticks(rotation=75)
plt.tight_layout()
plt.savefig('pipelineperf.png')
