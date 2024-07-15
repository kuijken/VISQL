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
	mkdir FITS BIAS CHARGE DARK FLAT SCIENCE GAIAPHOTOM
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

