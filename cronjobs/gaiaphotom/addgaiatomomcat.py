import numpy as np
from astropy.table import Table
from sys import argv

root=argv[1][:42]

momcat=Table.read(root+'_stars2_mom.cat',format='ascii')
frat  =Table.read(root+'_frat.txt',format='ascii')

z=momcat['XWORLD'] * 0. - 99.
momcat.add_columns([z,z],names=['G','BpRp'])
del z

for i in range(len(momcat)):
    infrat=frat['XWORLD'] == momcat['XWORLD'][i]
    if infrat.sum():
        ii=np.argmax(infrat)
        momcat['G'][i]=frat['G'][ii]
        momcat['BpRp'][i]=frat['BpRp'][ii]

momcat=momcat[momcat['G'] > 0.]

momcat.write(root+'_stars2_mom_col.cat',
             format='ascii.commented_header',
             overwrite=True)
