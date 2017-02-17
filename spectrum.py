# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 21:47:37 2017

@author: Alex
"""
import numpy as np

import matplotlib.pyplot as plt

def alpha(p):
    return 1/np.sqrt(1/p-1)

def lmbda(j,k,N,p):
    if j==0:
       return k/np.sqrt(N)+np.sqrt(N)*alpha(p)
    else:
        return (np.sin(np.pi*(k+1)*j/N)/np.sin(np.pi*j/N)-1)/np.sqrt(N)
    
def GD(q,k,N,p): #our spectrum
    s=0;
    for j in range (0,N-1):
        s=s+1/(q-lmbda(j,k,N,p))
    return s

def A(p,N):
    a = np.random.rand(N,N)  
    for i in range(N):
        for j in range(N):
            if a[i,j]<p:
                a[i,j]=(1/p-1)/np.sqrt(N)
            else:            
                a[i,j]=-1*alpha(p)/np.sqrt(N)
    return a;

def C(k,N):
    c = np.zeros((N,N)) 
    vec = np.zeros(N)
    
    for i in range(N):
        if ((1<=i<=k/2) or (i>=N-k/2)):
            vec[i]=1
    
    for i in range(N):
        for j in range(N):
            c[i,j]=vec[i-j]
    return c

def J(N):
    return np.ones((N,N))

def S(gamma,p,N,k):
    return A(p,N)+(gamma*C(k,N)+alpha(p)*J(N))/np.sqrt(N)

def plot(re,im):
    plt.figure(1)
    plt.subplot(211)
    plt.hist2d(re, im, bins=50)
    plt.colorbar()
    #plt.plot(values,'ro')
    plt.title('Eigenvalues of D')
    plt.xlabel('real(lambda)')
    plt.ylabel('imag(lambda)')
    #plt.axis([-1,N+1,-1,3])
    plt.show()
    #plt.show()
    
    plt.figure(2)
    plt.subplot(212)
    plt.plot(re, im, 'ro')
    return
    

k=6 #connectivity even
p=.000001 #0<p<1
N=1000 # size of our square matrix
gamma=1

eigvals=np.linalg.eigvals(S(gamma,p,N,k))
re = [z.real for z in eigvals]
im = [z.imag for z in eigvals]

plot(re,im)


