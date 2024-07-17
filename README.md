# VISQL

Quick-look analysis of Euclid-VIS raw fits files

K.Kuijken, 2023-2024


VIS quick-look raw fits data monitoring code
============================================

All code is wrapped by tcsh scripts.

Copy dot.p3start to .p3start in home directory to define a few p3 startup modules
The following python modules are imported by various scripts and need to be present:

<pre>
from astropy.coordinates import get_sun
from astropy.table import Table
from astropy.time import Time
from astroquery.vizier import Vizier
from datetime import date
from datetime import datetime
from glob import glob
from matplotlib import use
from os import environ
from os.path import basename
from sys import argv
import astropy.coordinates as coord
import astropy.io.ascii as asc
import astropy.io.fits as pf
import astropy.io.votable as vot
import astropy.units as u
import datetime
import datetime as dt
import json
import match
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import os 
import os, requests
import time
</pre>

the Gaiaphotom scripts use custom match.py and angsep.py from Henry Ferguson (included in the gaiaphotom folder).

---------

The code copies all fits files to a single directory, and soft-links the fits files to folders for the different types of data (bias, science, etc.)
All folders live under a single top-level directory, pointed to by the environment variable $viskom, which you need to configure by hand in the 'viskom' script.

Make sure that $viskom is set to the toplevel directory, and do
<pre>
	cd $viskom
	mkdir FITS BIAS CHARGE DARK FLAT SCIENCE TEST TEST/GAIAPHOTOM
</pre>

Then copy the contents of all cronjobs/* directories to the $viskom/src directory (do not make subdirectories, all code should live in $viskom/src/ ). 
Some files are duplicated, you should be able to select only the ones you want.

External executables that are needed to run these scripts include
<pre>
SExtractor, dfits, dos2unix, ds9 .
</pre>


The hierarchy of processing is as follows (so that e.g., catproc requires setyyyymm, dleas and sciencexy):

<pre>
setyyyymm, dleas
    sciencexy
        catproc
	    psf
                gaiaphotom
    darketc
    png

xmovie combines results from png, catproc and darketc.
</pre>
---------

The following jobs, which can be run as cronjobs, take care of various parts of the processing:

**cronjobdleas** -
     query the EAS for the LE1 fits files of the past 2 weeks, and download the FITS files not yet present

**cronjobsciencexy** -
     sort fits files by type (SCIENCE, DARK, CHARGE, BIAS, FLAT) based on header information,
     and make SExtractor catalogues for SCIENCE frames. There are 3 jobs that run concurrently:<br>
     science01 processes fits files whose yyyymmdd_hhmmss has ss beginning with 0 or 1 <br>
     science23 processes fits files whose yyyymmdd_hhmmss has ss beginning with 2 or 3 <br>
     science45 processes fits files whose yyyymmdd_hhmmss has ss beginning with 4 or 5 <br>

**cronjobcatproc** -
     process science SExtractor catalogues, selecting and measuring stars, and selecting/plotting cosmics

**cronjobdarketc** -
     process DARK, FLAT, BIAS fits files, a subset of what is done for SCIENCE frames

**cronjobpng** -
     make thumnails and zoom-in png images of all fits files; and run astrometry.net on science png's.

**cronjobpsf** -
     aggregate the star measurements from the science images (from catproc) and produce IQ scores

**cronjobxmovie** -
     correlate GOES solar Xray fluxes with OBSDATE, and plot CR distributions for SCIENCE and DARK.
     Then makes movies of the result, one showing proton counts and one showing GOES Xrays next to VIS.

**cronjobsetyyyymm** -
     sets a file '$scr/yyyymm' that contains the months to be processed.
     This job selects current month and last month if at most 3 days ago, but you can list any month(s)
     in this file if desired. The more months the more processing and checking will take place though.<br>
     For example if the file contains the lines
     <pre>
     202405
     202406
     </pre>
     then the files for those 2 months will be scanned and any missing catalogues, plots etc. will be generated.

**gaiaphotom** -
     match star catalogues to GAIA, make (VIS-G) vs (Bp-Rp) colour-colour diagram and deduce relative zeropoint
     aggregate results in various plots and tables (per exposure, per day, per quadrant, per ccd, etc)

---------

**Example output**

(see example_output folder)

For each SCIENCE frame:
<pre>
	In SCIENCE/ directory:
C_20240713_032840_W_52167482.VIS.bin_01_01.fits			- link to fits file
C_20240713_032840_W_52167482.VIS.bin_01_01.cat.gz		- SExtractor catalogue
C_20240713_032840_W_52167482.VIS.bin_01_01.png			- thumbnail image of full mosaic
C_20240713_032840_W_52167482.VIS.bin_01_01.wcs			- WCS for thumbnail from astrometry.net
C_20240713_032840_W_52167482.VIS.bin_01_01_cc.png		- 3x3 grid of full-resolution thumbnails
C_20240713_032840_W_52167482.VIS.bin_01_01_cl.png			[top,center,lower][left,center,right]
C_20240713_032840_W_52167482.VIS.bin_01_01_cr.png
C_20240713_032840_W_52167482.VIS.bin_01_01_lc.png
C_20240713_032840_W_52167482.VIS.bin_01_01_ll.png
C_20240713_032840_W_52167482.VIS.bin_01_01_lr.png
C_20240713_032840_W_52167482.VIS.bin_01_01_tc.png
C_20240713_032840_W_52167482.VIS.bin_01_01_tl.png
C_20240713_032840_W_52167482.VIS.bin_01_01_tr.png
C_20240713_032840_W_52167482.VIS.bin_01_01_cosmiX.png		- image of cosmic ray counts, around ADU for Xrays
C_20240713_032840_W_52167482.VIS.bin_01_01_cosmics.png		- ditto, wider ADU range
C_20240713_032840_W_52167482.VIS.bin_01_01_protons.png		- ditto, only high energies (protons)
C_20240713_032840_W_52167482.VIS.bin_01_01_rf.png		- flux v radius plot (for Xray movie)
C_20240713_032840_W_52167482.VIS.bin_01_01_chimney.jpg		- chimney plot
C_20240713_032840_W_52167482.VIS.bin_01_01_chimney2.png		- chimney plot with star sequende identified
C_20240713_032840_W_52167482.VIS.bin_01_01_chimney4.png		- chimney plot colour coded by ellipticity
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2.cat		- star catalogue (subset of full SExtractor catalogue)
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2.png		- zoomed images of PSF stars across the focal plane
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_mom.cat	- moment measurements of the stars in the star catalogue
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_mom.png	- plots of moments vs X and Y
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_mom_map3.png	- maps of R2, ellipticity, coma, trefoil, median-binned per CCD
C_20240713_032840_W_52167482.VIS.bin_01_01_iq.png		- mosaic image showing thumbnail, chimney, and ellipticity map
C_20240713_032840_W_52167482.VIS.bin_01_01_goes.png		- evolution of GOES X-ray flux in the hours around the exposure
	FROM GAIA MATCHING, IN TEST/GAIAPHOTOM/ directory:
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_gaia.cat		- catalogue of all Gaia stars within 0.5deg of nominam RA, DEC
C_20240713_032840_W_52167482.VIS.bin_01_01_x_gaia_dxdy.txt		- output from matching Gaia to star cat
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_x_gaia_dxdy.png	- plot of astrometry residuals vs Gaia
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_gaia.png		- plot of Gaia stars on outline of the FPA
C_20240713_032840_W_52167482.VIS.bin_01_01_frat.txt			- VIS/G flux ratio catalogue for all stars matched to Gaia
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_mom_col.cat		- star moments catalogue matched to Gaia, with G and Bp-Rp
C_20240713_032840_W_52167482.VIS.bin_01_01_stars2_fluxrat_bprp.png	- Plot of VIS-G vs Bp-Rp (for flux loss measurement)
</pre>

**Aggregated results on image quality:**

In the example_output/imagequality folder, for every month yyyymm (eg
202407), the following tables and plots with median IQ values per
exposure:
<pre>
IQ_timeline_202407.png                - plot of R2, e1, e2, tref1, tref2 vs. time (median of indiv. exposures is plotted)
psfmoments_202407.png                 - plot of combinations of e,R2,coma,trefoil
psfmoments_202407.txt                 - Table of median IQ parameters per exposure (redmine wiki format)
psfmoments_202407_le1_score.csv       - Table of R2, e1, e2 and IQ score, with EAS name of exposure (csv format)
psfmoments_202407_le1_score.txt       - Table of R2, e1, e2 and IQ score, with EAS name of exposure (wiki format)
psfmoments_202407_score.csv           - Table of R2, e1, e2 and IQ score (csv format)
psfmoments_202407_score.txt           - Table of R2, e1, e2, coma1,2, trefoil1,2 (wiki format)
</pre>

**Aggregated results on Gaia photometry comparison:**

In the example_output/gaiaphotometry folder, a number of aggregated
plots and tables of the average fitted VIS/Gaia flux ratio at three
Bp-Rp colours, per day and over longer periods. Most have 'long' and
'short' versions, for normal science (560) and short science (89)
second exposures.

<pre>
EVERY DAY:
iceplot_FPA_240712_560.52.png             - VIS/Gaia at Bp-Rp=1.75 hexbinned over FPA, data for one day 
vis-over-gaia-20240711-ccd_short.txt      - table of day-averaged mean VIS/Gaia vs Bp-Rp fits, binned by quadrant, ccd, and 2x2 CCD block 

EVERY WEEK:
iceplot_FPA_wk240711-240718_560.52.png    - VIS/Gaia at Bp-Rp=1.75 hexbinned over FPA, data for one week 

EVERY MONTH:
VISGaia_FPA_SHORT_202407.png              - VIS/Gaia vs Bp-Rp binned per CCD, data for entire month
iceplot_FPA_202407_560.52.png             - VIS/Gaia at Bp-Rp=1.75 hexbinned over FPA, data for entire month
vis-over-gaia-202407-3x3_short.txt        - table of 3x3 binned mean VIS/Gaia vs Bp-Rp fits, for all exposures in the month 
vis-over-gaia-202407_long.txt             - table of FPA-averaged mean VIS/Gaia vs Bp-Rp fits, for all exposures in the month 
visgaiacolour_LONG_202407.txt             - table of month-averaged per-quadrant and per-ccd mean VIS/Gaia vs Bp-Rp fits 

LONGER TIMELINES:
iceplot_diag_short.png                    - time evolution of flux loss on six diagonal CCDs, daily averages, from March 2024
iceplot_xdiag_short.png                   - time evolution of flux loss on six cross-diagonal CCDs, daily averages, from March 2024
iceplot_short.png                         - time evolution of mean flux loss, averaged per exposure, three Bp-Rp colours, full timeline
iceplot_short2.png                        - time evolution of mean flux loss, averaged per exposure, three Bp-Rp colours, from March 2024
iceplot_short3.png                        - time evolution of mean flux loss, averaged per exposure, three Bp-Rp colours, from Jun 2024
iceplot_selfcal_long.png                  - as above, selfcal field data only
iceplot_selfcal_long2.png                 - as above, selfcal field data only
iceplot_selfcal_long3.png                 - as above, selfcal field data only
</pre>
