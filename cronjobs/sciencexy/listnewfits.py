from glob import glob
from os import environ,rename
from sys import argv

# check for 'clean' command to remove corrupted lines from fitslist.txt
writecleanfitslist=False
if len(argv)>1:
    if argv[1]=='clean':
        writecleanfitslist=True

# read existing fitslist.txt file
fitslist = [x for x in open('fitslist.txt','r')]

# print out any illegal lines
print('Checking for corrupted lines in fitslist.txt:')
badlines=[]
for x in fitslist:
    xs=x.split()
    if (xs[0]!='#') & ((len(xs)!=5) | (len(xs[0])!=47)):
        print(x,end='')
        badlines.append(x)
        
# internally only keep the legal lines
print(len(fitslist),'lines in fitslist.txt were found.')
for x in badlines:
    fitslist.remove(x)
print(len(fitslist),'lines in fitslist.txt were OK.')
# make a set of the fits file names in the list
fitsnamesinlist=set([x.split()[0] for x in fitslist])

# if requested, write out a clean version of fitslist
if writecleanfitslist:
    rename('fitslist.txt','fitslist.txt.SAVED')
    f=open('fitslist.txt','w')
    for x in fitslist:
        print(x,file=f,end='')
    f.close()
    print('Clean version of fitslist written out and old one moved to fitslist.txt.SAVED')
    

# make list of the fits files in local directory for the months in $scr/yyyymm
yyyymm=[x[:-1] for x in open(environ['viskom']+'/scr/yyyymm','r')]
fitsfiles=[]
for ym in yyyymm:
    fitsfiles += glob('C_'+ym+'*.fits')
fitsfiles.sort()

# make a list of the fits files that were not yet in the list
newfitsfiles=[x for x in fitsfiles if x not in fitsnamesinlist]
newfitsfiles.sort()
# and write out to tmpnewfitsfiles.txt
f=open('tmpnewfitsfiles.txt','w')
for x in newfitsfiles:
    print(x,file=f)
f.close()
    


