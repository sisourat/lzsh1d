import sys 
import numpy as np
import mydistrib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

finp = open(sys.argv[1],"r")
emin = float(sys.argv[2])
emax = float(sys.argv[3])
de = float(sys.argv[4])
r = []
for l in finp:
    w = l.split()
    r.append(float(w[0]))
    r.append(float(w[1]))

y, x = np.histogram(r, bins=np.arange(emin,emax,de))
ntot = len(r)
for i in range(len(y)):
    if(y[i]>0):
        sig = float(y[i])/ntot
        dsig = sig*np.sqrt( float(ntot-y[i]) / (ntot*y[i]))
        print  x[i], y[i]#,sig, dsig
    else:
        print  x[i], 0.0, 0.0










