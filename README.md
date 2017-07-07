# Turing-Patterns
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:53:15 2017

@author: Sakshi Khaiwal
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 18:55:04 2017

@author: Sakshi Khaiwal
"""

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as patches
from matplotlib.widgets import Slider
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy as dc

L=2
DA=0.00001  #Diffusion constant for A
DB=0.001    #Diffusion constant for A
K=80      #Number of compartments
h=float(L)/K      #length of each compartment
da=DA/(h*h) #rate of diffusion of A
db=DB/(h*h) #rate of diffusion of B 
n_reactions= 100000 #number of reactions
k1=0.000001  #Degradation rate of A
k2=1     #rate of production of moleculesA
k3=0.002      #Degradation rate of A
k4=3      #rate of production of moleculesA
A=np.zeros(K) #number of molecules of A in each compartment
B=np.zeros(K) #number of molecules of B in each compartment

x=np.zeros(K)   #postions
a1 = np.zeros(K) #propensity function for diffusion of A till K-1
a2 = np.zeros(K) #propensity function for diffusion of A till K
b1 = np.zeros(K) #propensity function for diffusion of B till K-1
b2 = np.zeros(K) #propensity function for diffusion of B till K
a3 = np.zeros(K) #propensity function for degradation of A 
a4 = np.zeros(K) #propensity function for production of A
a5 = np.zeros(K) #propensity function for degradation of B
a6 = np.zeros(K) #propensity function for production of B
Alis = np.zeros((n_reactions//1000+1,K)) 
Blis = np.zeros((n_reactions//1000+1,K))
t=np.zeros(n_reactions)
time=0
for m in range(K):  #initial condition
    A[m]=100
    B[m]=40
     

Alis[0,:]=A[:]
Blis[0,:]=B[:]
c1=0

tlis = [0,]
f=open("turingpattern.txt","+w")
for i in range (1,n_reactions):
    r1=np.random.uniform(0,1)
    r2=np.random.uniform(0,1)
    a1[0]=0
    b1[0]=0
    for k in range(1,K):
        a1[k]= a1[k-1]+A[k]*da
        b1[k]= b1[k-1]+B[k]*db
    
    a2[-1]=0
    b2[-1]=0
    for k in range (K-1):
        a2[k]= a2[k-1]+A[k]*da
        b2[k]= b2[k-1]+B[k]*db
          

    a3[-1]=0
    a4[-1]=0
    a5[-1]=0
    a6[-1]=0
    for k in range(K):
        x[k] = k*h
        a3[k]=0.5* A[k]*k1*(A[k]-1)*B[k]+a3[k-1]
        a4[k]= k2 + a4[k-1]
        a5[k]= A[k]*k3 + a5[k-1]
        a6[k]= k4 + a6[k-1]  
    
    a1f = a1[K-1]
    a2f = a2[K-2]
    b1f = b1[K-1]
    b2f = b2[K-2]
    a3f = a3[K-1]
    a4f = a4[K-1]
    a5f = a5[K-1]
    a6f = a6[K-1]
    #if i==100: print(a1f,a2f,a3f,a4f,a5f,a6f)
    a0 = a1f+a2f+a3f+a4f+a5f+a6f+b1f+b2f
    tau = (1/a0)*(np.log(1/r1))
    time = time +tau
    #t[i]=time
    #A[:,i]=A[:,i-1]
    #B[:,i]=B[:,i-1]
       #break
       #print(time)
    if (0<r2<(a2f/a0)):
        for j in range (K-1):
            if (r2< (a2[j]/a0)):
                A[j]=A[j]-1
                A[j+1]=A[j+1]+1
                break
    elif(r2<(a2f+a1f)/a0):
        for j in range (1,K):
            if (r2< ((a2f+a1[j])/a0)):
                A[j]=A[j]-1
                A[j-1]=A[j-1]+1
                break
    elif (r2<((a2f+a1f+b2f)/a0)):
        for j in range (K-1):
            if (r2< ((a1f+a2f+b2[j])/a0)):
                B[j]=B[j]-1
                B[j+1]=B[j+1]+1
                break
    elif(r2<(a2f+a1f+b1f+b2f)/a0):
        for j in range (1,K):
            if (r2< ((a1f+b2f+a2f+b1[j])/a0)):
                B[j]=B[j]-1
                B[j-1]=B[j-1]+1
                break
    elif(r2<(a2f+a1f+b1f+b2f+a3f)/a0):
        for j in range (K):
            if (r2< ((a2f+a1f+b1f+b2f+a3[j])/a0)):
                A[j]=A[j]+1
                B[j]=B[j]-1
                
                break
    elif(r2<(a2f+a1f+b1f+b2f+a3f+a4f)/a0):
        c1+=1
        for j in range (K):
            if (r2< ((a2f+a1f+b1f+b2f+a3f+a4[j])/a0)):
                A[j]=A[j]+1
                break
    elif(r2<(a2f+a1f+b1f+b2f+a3f+a4f+a5f)/a0):
        for j in range (K):
            if (r2< ((a2f+a1f+b1f+b2f+a3f+a4f+a5[j])/a0)):
                A[j]=A[j]-1
                break
    elif(r2<(a2f+a1f+b1f+b2f+a3f+a4f+a5f+a6f)/a0):
        for j in range (K):
            if (r2< ((a2f+a1f+b1f+b2f+a3f+a5f+a6[j])/a0)):
                B[j]=B[j]+1
                break
    """if i%1000==0:
        #print(A)
        tlis.append(time)
        #print(i//1000)
        Alis[i//1000,:]=A[:]
        Blis[i//1000,:]=B[:]#"""
print(time)    
for j in range(K):
        f.write("%d\r\n"%(A[j])+"  "+"%d\r\n"%(B[j]))
f.close()
"""   
print(c1)
Alis=np.array(Alis)
Blis=np.array(Blis)
tlis=np.array(tlis)
#"""
"""
fig, ax = plt.subplots(1, figsize=(8,6))

ax.set_xlim(x.min(), x.max())
#ax.set_ylim(A.min(), A.max())
ax.set_ylim(min(Alis.min(),Blis.min()), max(Alis.max(),Blis.max()))
aplot, = ax.step(x, Alis[0], where='mid')
bplot, = ax.step(x, Blis[0], where='mid')
ax.set_xlabel('cell positions')
ax.set_ylabel('number of A molecules')

divider = make_axes_locatable(ax)
sax = divider.append_axes("top",size="10%",pad=0.2,)
sl = Slider(sax, 't', tlis.min(), tlis.max(), valinit=tlis.min())

# update function
def update(val):
    ti = (abs(t-val)).argmin()
    aplot.set_ydata(Alis[ti])
    bplot.set_ydata(Blis[ti])
    plt.draw()
sl.on_changed(update)

plt.show()
"""
plt.bar(x,A,0.025)
plt.show()
#"""
