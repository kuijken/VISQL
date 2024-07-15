from sys import argv
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use

use('agg')

fitsfile = argv[1]
catfile  = fitsfile[:-5]+"_stars2.cat"
aspar    = 10     # asinh parameter - linear up to this level, then log
if len(argv)>2: aspar=float(argv[2])

x,y,flx,flxerr,r,a,b,pa,id,flg,rim = np.loadtxt(catfile,unpack=True)
keep = (b>1.1)
x=x[keep]
y=y[keep]
flx=flx[keep]
flxerr=flxerr[keep]
r=r[keep]
a=a[keep]
b=b[keep]
pa=pa[keep]
id=id[keep]
flg=flg[keep]
rim=rim[keep]
starid=np.array([i for i in range(keep.sum())])

#read in all extensions CRPIX1,2, NAXIS 1,2, en CD1_1, CD2_2 values
nout=400
cout=nout//2
out=np.zeros((nout,nout))
degwide=0.85
plt.figure(figsize=(0.04*nout,0.04*nout),dpi=100)
w=10 # thumnail extent from center pixel
f=pf.open(fitsfile)
for ext in f[1:]:
    extname=ext.header['EXTNAME']
    crpix1=ext.header['CRPIX1']
    crpix2=ext.header['CRPIX2']
    naxis1=ext.header['NAXIS1']
    naxis2=ext.header['NAXIS2']
    cd1_1 =ext.header['CD1_1']
    cd2_2 =ext.header['CD2_2']
    crval1=ext.header['CRVAL1']
    crval2=ext.header['CRVAL2']
    #find all stars on a given extension
    i=(x-crval1)/cd1_1+crpix1
    j=(y-crval2)/cd2_2+crpix2
    # take a margin of 100 pix from the edges
    onchip = (i>100) & (i<naxis1-100) & (j>100) & (j<naxis2-100)
    if onchip.sum()>0:
#        print(onchip.sum(),'stars found on extension',extname)
        pixels=f[extname].data
        bg=np.median(pixels[100:-100,100:-100])
        for star in starid[onchip]:
            istar,jstar=int(i[star]+0.5),int(j[star]+0.5)
            zoom=pixels[jstar-w:jstar+w-1,istar-w:istar+w-1]-bg
#            plt.imshow(np.arcsinh(zoom/aspar)*aspar,vmin=0,origin='lower')
#            plt.show()
            iout=int(cout+x[star]*nout/degwide)
            jout=int(cout+y[star]*nout/degwide)
#            only paint in the star if it does not overlap with an earlier one.
            if out[jout-w:jout+w-1,iout-w:iout+w-1].sum()==0:
                out[jout-w:jout+w-1,iout-w:iout+w-1]+=zoom
plt.imshow(np.arcsinh(out/aspar)*aspar,vmin=0.,origin='lower',
           interpolation='none')
plt.colorbar(label='%0.f arcsinh(ADU/%.0f)' % (aspar,aspar))
plt.title(fitsfile)
plt.savefig(fitsfile[:-5]+'_stars2.png')
#plt.show()
