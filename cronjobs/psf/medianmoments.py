import numpy as np
from sys import argv

print('#|   catalogue                                                 |'+
          '      R2  |    COMA1 |  COMA 2 |    ELL 1 |   ELL 2 |'+
          '   TREF 1 |  TREF 2 |')
for momcat in argv[1:]:
    r2,e1,e2,t1,t2,c1,c2=np.loadtxt(momcat,unpack=True,
                                        usecols=[12,13,14,15,16,17,18])
    print('| %-60s | %8.3f | %8.4f |%8.4f | %8.4f |%8.4f | %8.4f |%8.4f |'
              % tuple([momcat]+[np.median(x) for x in [r2,c1,c2,e1,e2,t1,t2]]),
              )

