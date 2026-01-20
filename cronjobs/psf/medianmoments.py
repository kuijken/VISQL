import numpy as np
from sys import argv
from astropy.table import Table

print('#|   catalogue                                                 |'+
          '      R2  |    COMA1 |   COMA 2 |    ELL 1 |    ELL 2 |'+
          '   TREF 1 |   TREF 2 |')

for momcat in argv[1:]:
    try:
        t=Table.read(momcat,format='ascii')
        r2=t['R2']
        e1=t['e1']
        e2=t['e2']
        t1=t['tr1']
        t2=t['tr2']
        c1=t['com1']
        c2=t['com2']
        if len(r2)==0:
            r2,e1,e2,t1,t2,c1,c2=[9.9], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999]
    except:
        r2,e1,e2,t1,t2,c1,c2=[9.9], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999]

    print('| %-60s | %8.3f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f |'
              % tuple([momcat]+[np.median(x) for x in [r2,c1,c2,e1,e2,t1,t2]]),
              )








# for momcat in argv[1:]:
#     try:
#         r2,e1,e2,t1,t2,c1,c2=np.loadtxt(momcat,unpack=True,
#                                         usecols=[12,13,14,15,16,17,18])
#     except:
#         r2,e1,e2,t1,t2,c1,c2=[99], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999], [0.9999]

