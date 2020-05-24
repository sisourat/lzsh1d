import sys
import numpy as np

r1 = np.linspace(2.0, 202.0, num=40)
r2 = np.linspace(2.0, 202.0, num=40)
th = np.linspace(0.0, 90.0, num=3)

print 40, 40, 3
for i1 in r1:
  for i2 in r2:
    for it in th:
      print i1,i2,it,1.0/i1**4+1.0/i2**2
