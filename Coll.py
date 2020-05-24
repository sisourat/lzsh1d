import sys
from shutil import copyfile
import os
import numpy as np
import time

fcurves=sys.argv[1]
rf=sys.argv[2]
ista=sys.argv[3]

ecoll=10.0 # eV
ecollau=ecoll/27.211

m=1.0
z0=-40.0
bmin=2.0    # au
bmax = 38.0 # au
nb = 72
blist=np.linspace(bmin,bmax,num=nb)
ntraj = 200
v=np.sqrt(2.0*ecollau/(m*1836.15))

ib=0
bpb1=np.zeros(nb+1)
bpb2=np.zeros(nb+1)
bpb3=np.zeros(nb+1)
for b in blist:
 wsave = []
 ib+=1
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
  os.system('/home/nico/Workspace/Propag_LZSH_1D/dyn < traj > a')
#  os.system('/home/nico/Results/Patrik/lzsh_sources/dyn < traj > a')

#  cmd='cat a' 
#  os.system(cmd)
#  cmd='cp a traj'+str(j) 
#  os.system(cmd)
#  cmd='cp fort.* traj'+str(j) 
#  os.system(cmd)

  fb=open("a","r")
  l=fb.readline()
  w=l.split()
  if(w[0]=='error'):
   print b
   print l
   sys.exit()
   
#  print l
  wsave.append(w[0])
  if(not(w[0]==ista)):
   pb+=1.0
  if(w[0]=='1'):
   pb1+=1.0
  if(w[0]=='2' or w[0]=='3'):
   pb2+=1.0
  if(w[0]=='4' or w[0]=='5' or w[0]=='6' or w[0]=='7' or w[0]=='9'):
   pb3+=1.0
 print b,0.28*2.0*np.pi*b*pb/float(ntraj),0.28*2.0*np.pi*b*pb1/float(ntraj),0.28*2.0*np.pi*b*pb2/float(ntraj),0.28*2.0*np.pi*b*pb3/float(ntraj)#,pb,pb/float(ntraj)
 bpb1[ib]=2.0*np.pi*b*pb1/float(ntraj)
 bpb2[ib]=2.0*np.pi*b*pb2/float(ntraj)
 bpb3[ib]=2.0*np.pi*b*pb3/float(ntraj)
 fout=open("out_b_"+str(b)+".txt","w")
 print >> fout, b,pb1/float(ntraj),pb2/float(ntraj), pb3/float(ntraj)
 print >> fout, len(wsave)
 for ww in wsave:
  print >> fout, ww
