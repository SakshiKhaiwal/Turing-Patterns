# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 15:53:28 2017

@author: Sakshi Khaiwal
"""
import numpy as np
import matplotlib.pyplot as plt

DA=0.00001  #Diffusion constant
DB=0.001
finaltime=1000# four minutes 

L = 1            # Total length
K=40             # number of compartments
h=float(L)/K     # Length of each compartment
da=DA/(h*h)
db=DB/(h*h)
x = np.zeros(K)
A = np.zeros(K)
B = np.zeros(K)
a1 = np.zeros(K)
b1 = np.zeros(K)
a2 = np.zeros(K)
b2 = np.zeros(K)

'''for i in range(K):
    A[i]=200
    B[i]=75'''
    
A[20]=500
A[21]=500
B[20]=0
time=0

while (time<finaltime):
       r1=np.random.uniform(0,1)
       r2=np.random.uniform(0,1)
       a1 = np.zeros(K)
       b1 = np.zeros(K)
       a2 = np.zeros(K)
       b2 = np.zeros(K)
       for i in range (K-1):
           x[i] = i*h
           a2[i]= a2[i-1]+A[i]*da
           b2[i]= b2[i-1]+B[i]*db
       for i in range(1,K):
           a1[i]= a1[i-1]+A[i]*da
           b1[i]= b1[i-1]+B[i]*db
       x[K-1] = h*(K-1)
       a1f = a1[K-1]
       b1f = b1[K-1]
       a2f = a2[K-2]
       b2f = b2[K-2]
       a0 = a1f+a2f+b1f+b2f
       tau = (1/a0)*(np.log(1/r1))
       time = time +tau
       #break
       #print(time)
       if (0<r2<(a2f/a0)):
           #print 1,
           for j in range (K-1):
               if (r2< (a2[j]/a0)):
                   #print j,
                   A[j]=A[j]-1
                   A[j+1]=A[j+1]+1
                   break
       elif(r2<(a2f+a1f)/a0):
           #print 2,
           for j in range (1,K):
               if (r2< ((a2f+a1[j])/a0)):
                   #print j
                   A[j]=A[j]-1
                   A[j-1]=A[j-1]+1
                   break
       elif (r2<(a1f+a2f+b2f)/a0):
           #print 3,
           for j in range (K-1):
               if (r2< (a1f+a2f+b2[j])/a0):
                   B[j]=B[j]-1
                   B[j+1]=B[j+1]+1
                   break
       else:
           #print 4,
           for j in range (1,K):
               if (r2< (a1f+a2f+b2f+b1[j])/a0):
                   B[j]=B[j]-1
                   B[j-1]=B[j-1]+1
                   break
       print (time)
plt.figure() 
#plt.subplot(2,1,1)                   
plt.bar(x,A,0.025)
'''
plt.subplot(2,1,2)                   
plt.plot(x,B)'''

plt.show()                          