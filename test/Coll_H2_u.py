import sys
from shutil import copyfile
import os
import numpy as np
import time


# input in command line
fcurves=sys.argv[1]
rf=sys.argv[2]
ista=sys.argv[3]

# hard-coded parameters (should be adjusted for each system)

ecoll=1.0 # eV
ecollau=ecoll/27.211
m=1.0
v=np.sqrt(2.0*ecollau/(m*1836.15))

z0=-40.0
bmin=2.0    # au
bmax = 38.0 # au
nb = 72
ntraj = 200
blist=np.linspace(bmin,bmax,num=nb)
kerbins = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0]
dker=0.2


# initialize some arrays
ib=0
bpbtot=np.zeros(nb)
bpb1=np.zeros(nb)
bpb2=np.zeros(nb)
bpb3=np.zeros(nb)

stot = 0.0
s1 = 0.0
s2 = 0.0
s3 = 0.0

skertot = np.zeros((nb,len(kerbins)-1))
sker1 = np.zeros((nb,len(kerbins)-1))
sker2 = np.zeros((nb,len(kerbins)-1))
sker3 = np.zeros((nb,len(kerbins)-1))

# dynamics

for b in blist:
 wsave = []
 ker1 = []
 ker2 = []
 ker3 = []
 kertot = []
 pb = 0.0
 pb1 = 0.0
 pb2 = 0.0
 pb3 = 0.0
 for j in range(ntraj):
  finp=open("traj","w")
  print >> finp, fcurves,rf,ista
  print >> finp, m,0.0,0.0,0.0,0.0,0.0,0.0
  print >> finp, m,0.0,b,z0,0.0,0.0,v
  finp.close()
  os.system('/home/nico/Workspace/Propag_LZSH_1D/dyn_stable < traj > a')

#  os.system('/home/nico/Results/Patrik/lzsh_sources/dyn < traj > a')
#  cmd='cp a traj'+str(j) 
#  os.system(cmd)
#  cmd='cp fort.* traj'+str(j) 
#  os.system(cmd)

# analyze each traj.
  fb=open("a","r")
  l=fb.readline()
  w=l.split()
  if(w[0]=='error'):
   print b
   print l
   sys.exit()
   
  wsave.append(w[0])
  if(not(w[0]==ista)):
   pb+=1.0
   kertot.append(float(w[1])-ecoll)
  if(w[0]=='1' or w[0]=='2'):
   pb1+=1.0
   ker1.append(float(w[1])-ecoll)
  if(w[0]=='3' or w[0]=='4' or  w[0]=='5' or w[0]=='6'):
   pb2+=1.0
   ker2.append(float(w[1])-ecoll)
  if(w[0]=='8'):
   pb3+=1.0
   ker3.append(float(w[1])-ecoll)

# analysis for each b
 bpbtot[ib]=2.0*np.pi*b*pb/float(ntraj)*0.28
 bpb1[ib]=2.0*np.pi*b*pb1/float(ntraj)*0.28
 bpb2[ib]=2.0*np.pi*b*pb2/float(ntraj)*0.28
 bpb3[ib]=2.0*np.pi*b*pb3/float(ntraj)*0.28

 kertot_bin, kertot_bin_edges = np.histogram(kertot,bins=kerbins,density=False)
 skertot[ib] = 2.0*np.pi*b*0.28*kertot_bin/float(ntraj)/dker

 if(not(len(ker1)==0)):
   ker1_bin, ker_bin_edges = np.histogram(ker1,bins=kerbins,density=False)
 else:
   ker1_bin = np.zeros(len(kerbins)-1)
 sker1[ib] = 2.0*np.pi*b*0.28*ker1_bin/float(ntraj)/dker

 if(not(len(ker2)==0)):
   ker2_bin, ker_bin_edges = np.histogram(ker2,bins=kerbins,density=False)
 else:
   ker2_bin = np.zeros(len(kerbins)-1)
 sker2[ib] = 2.0*np.pi*b*0.28*ker2_bin/float(ntraj)/dker

 if(not(len(ker3)==0)):
   ker3_bin, ker_bin_edges = np.histogram(ker3,bins=kerbins,density=False)
 else:
   ker3_bin = np.zeros(len(kerbins)-1)
 sker3[ib] = 2.0*np.pi*b*0.28*ker3_bin/float(ntraj)/dker

 fout=open("out_H2_u_b_"+str(b)+"_ecoll_"+str(ecoll)+"eV_ntraj_"+str(ntraj)+".txt","w")

 print b,0.28*2.0*np.pi*b*pb/float(ntraj),0.28*2.0*np.pi*b*pb1/float(ntraj),0.28*2.0*np.pi*b*pb2/float(ntraj),0.28*2.0*np.pi*b*pb3/float(ntraj)
 print >> fout, b,pb1/float(ntraj),pb2/float(ntraj), pb3/float(ntraj)
 print >> fout, len(wsave)

 for k in range(len(kertot_bin)):
     print >> fout,(kertot_bin_edges[k]+kertot_bin_edges[k+1])/2.0, kertot_bin[k]*0.28*2.0*np.pi*b,ker1_bin[k]*0.28*2.0*np.pi*b, ker2_bin[k]*0.28*2.0*np.pi*b, ker3_bin[k]*0.28*2.0*np.pi*b

# for ww in wsave:
#  print >> fout, ww
 ib+=1

# total cross sections and ker spectra (integrated over b)
print
print np.trapz(bpbtot,x=blist),np.trapz(bpb1,x=blist),np.trapz(bpb2,x=blist),np.trapz(bpb3,x=blist)
print
fker=open("ker_H2_u_ecoll_"+str(ecoll)+"eV_ntraj_"+str(ntraj)+".txt","w")
for k in range(len(kertot_bin)):
    sker = skertot[:,k]
    print (kertot_bin_edges[k]+kertot_bin_edges[k+1])/2.0,np.trapz(sker,x=blist),
    print >> fker, (kertot_bin_edges[k]+kertot_bin_edges[k+1])/2.0,np.trapz(sker,x=blist),
    sker = sker1[:,k]
    print np.trapz(sker,x=blist),
    print >> fker, np.trapz(sker,x=blist),
    sker = sker2[:,k]
    print np.trapz(sker,x=blist),
    print >> fker, np.trapz(sker,x=blist),
    sker = sker3[:,k]
    print np.trapz(sker,x=blist)
    print >> fker, np.trapz(sker,x=blist)
#    for ib in range(len(blist)):
#       print blist[ib],(kertot_bin_edges[k]+kertot_bin_edges[k+1])/2.0,sker[ib]
#    print

    
