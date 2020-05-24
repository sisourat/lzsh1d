import sys
import numpy as np

e = np.array([-1.1034263742, -0.5493228411, -0.4766498111, -0.4761058109, -0.4754815269, -0.4517066386, -0.4512595793, -0.4508872360, -0.4403046702])

rarr = np.linspace(0.005,0.45,12)

for r in rarr:
 print r, ' '.join(map(str,e+1.0/np.power(r,0.5)-1/np.power(0.5,0.5)))
