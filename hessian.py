#
#  just type in command line : python hessian.py yourfile.dat
#  and it will print the frequencies on the screen
#  and the normal modes in normalcoord.txt
#

import sys
import numpy as np
from numpy import linalg as LA
import mydistrib

#######################################################
def convert(data):
         tempDATA = []
         for i in data:
                 tempDATA.append(float(i)/0.529)
         return np.array(tempDATA)
#######################################################
class atomicdata:
	def _init_(self):
		self.at = 'A'
		self.charge = 0.0
		self.mass = 0.0
		self.xyz = [0.0,0.0,0.0]

	def initMass(self):
		if(self.at=='C'):
			 self.mass = 12.0000*1836.15
		elif(self.at=='N'):
			 self.mass = 14.00307*1836.15
		elif(self.at=='O'):
			 self.mass = 15.99491*1836.15
		elif(self.at=='LI'):
			 self.mass = 7.0160*1836.15
		elif(self.at=='CL'):
			 self.mass = 34.96885*1836.15
		elif(self.at=='H' or self.at=='HYDROGEN'):
			 self.mass = 1.00782*1836.15
		elif(self.at=='Si' or self.at=='SI'):
			 self.mass = 27.97693*1836.15
		elif(self.at=='S' or self.at=='SI'):
			 self.mass = 32.066*1836.15
		else:
			print 'initMass not defined'
			sys.exit()
#######################################################

# reads FORCE CONSTANT matrix and Equilibrium Geometry from file.dat

ntraj = 2000

f = file(sys.argv[1],"r")

fQ = file('normalcoor.txt',"w")

w = ['0','0','0','0']
force = []
atoms = []
atom = atomicdata()
stop = 'n'

print 'READS DATA ON THE GS'

while(stop=='n'):
	line = f.readline()
	w = line.split()
	if(len(w)>0 and w[0]=='$HESS'):
		line = f.readline()
		while(not(line[0:5]==' $END')):
			line = f.readline()
			if(not(line[0:5]==' $END')):
				for i in range((len(line)-6)/15):
					j=i*15+5
					k=j+15
					force.append(float(line[j:k]))
#					print float(line[j:k])
		stop='y'
		print 'FULL FORCE MATRIX READ'
	if(len(w)>0 and w[0]=='$DATA'):
		while(not(w[0]=='$END')):
			line = f.readline()
			w = line.split()
			if(not(w[0]=='$END')):
				atom = atomicdata()
				atom.at = w[0]
				atom.charge = w[1]
				atom.xyz = convert(w[2:5])
				atom.initMass()
				atoms.append(atom)
		print 'ATOMIC DATA READ'

# builds the HESSIAN matrix

print 'BUILDS THE HESSIAN',len(force)

hess = np.zeros((3*len(atoms),3*len(atoms)))
for i in range(3*len(atoms)):
	for j in range(3*len(atoms)):
		ik = j + i*3*len(atoms)
		ia = i/3
		ja = j/3
		hess[i,j] = force[ik]/np.sqrt(atoms[ia].mass*atoms[ja].mass)

nmodes = 3*len(atoms)-6

# diagonalizes the hessian matrix
freq, v = LA.eig(hess)

print 'FREQUENCIES (CM**-1)'
print np.sqrt(freq[0:nmodes])*219474
print np.sqrt(freq[0:nmodes])
omega = np.sqrt(freq[0:nmodes])

p = v
pinv = LA.inv(p)

# prints the normal modes in normalcoord.txt
vec = np.transpose(v)
for i in range(nmodes):
        vv = vec[i]
        print >> fQ, np.sqrt(freq[i])*219474
	for j in range(len(atoms)):	
		print >> fQ, atoms[j].at,vec[i][j*3:j*3+3]

# transfrom back to cartesian coord. 

cr = []
cp = []
for i in range(nmodes):
    cr.append(np.sqrt(1.0/(2.0*omega[i])))
    cp.append(np.sqrt(omega[i]))
    print i,np.sqrt(1.0/(2.0*omega[i])), np.sqrt(omega[i])

Q0 = np.zeros(3*len(atoms))
PQ0 = np.zeros(3*len(atoms))
for itraj in range(ntraj):
     ftraj = open('traj_'+str(itraj)+'.init','w')
     for i in range(nmodes):
        Q0[i] = mydistrib.rand_normal(0.0,cr)
        PQ0[i] = mydistrib.rand_normal(0.0,cp)
#     print PQ0[0]
     chi=np.dot(np.transpose(pinv),Q0)
     pchi=np.dot(np.transpose(pinv),PQ0)
#     xyz = []
     for j in range(len(atoms)):
	k=3*j
        xa = atoms[j].xyz[0]+chi[k]/np.sqrt(atoms[j].mass)
        ya = atoms[j].xyz[1]+chi[k+1]/np.sqrt(atoms[j].mass)
        za = atoms[j].xyz[2]+chi[k+2]/np.sqrt(atoms[j].mass)
        pxa = pchi[k]#/np.sqrt(atoms[j].mass)
        pya = pchi[k+1]#/np.sqrt(atoms[j].mass)
        pza = pchi[k+2]#/np.sqrt(atoms[j].mass)
#        xyz.append([xa,ya,za])
	print >> ftraj, atoms[j].mass, xa, ya, za, pxa/atoms[j].mass, pya/atoms[j].mass, pza/atoms[j].mass, atoms[j].at #, (pxa**2+pya**2+pza**2)#*atoms[j].mass
     ftraj.close()
