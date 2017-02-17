# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 21:47:37 2017

@author: Alex
"""
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

def alpha(p):
    return 1/np.sqrt(1/p-1)

def lmbdaD(j,k,N,p):
    if j==0:
       return k/np.sqrt(N)+np.sqrt(N)*alpha(p)
    else:
        return (np.sin(np.pi*(k+1)*j/N)/np.sin(np.pi*j/N)-1)/np.sqrt(N)
    
def GD(q,k,N,p): #our spectrum
    s=0;
    for j in range (0,N-1):
        s=s+1/(q-lmbdaD(j,k,N,p))
    return s

def A(p,N):
    a = np.random.rand(N,N)  
    for i in range(N):
        for j in range(N):
            if a[i,j]<p:
                a[i,j]=(1/p-1)/np.sqrt(N)
            else:            
                a[i,j]=-1*alpha(p)/np.sqrt(N)
    #print(a)
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

def D(gamma,p,k,N):
    d= (gamma*C(k,N)+alpha(p)*J(N))/np.sqrt(N)
    #print(d)
    return d
    
def S(gamma,p,N,k):
    s=A(p,N)+D(gamma,p,k,N)
   # print(s)
    return s

def plot(re,im,b,N,p,k,mat):
    plt.figure(1)
    plt.subplot(211)
    plt.hist2d(re, im, bins=b,norm=LogNorm())
    plt.colorbar()
    title='Eigenvalues of {3},N={0},p={1},k={2}'.format(N,p,k,mat)
    filename='hist({0}=N)({1}=p)({2}=k)({3}=mat).png'.format(N,p,k,mat)
    plt.title(title)
    plt.xlabel('real(lambda)')
    plt.ylabel('imag(lambda)')
    plt.savefig(filename, bbox_inches='tight')
    
    plt.figure(2)
    plt.subplot(212)
    plt.plot(re, im,'bo',markersize=5)
    plt.title(title)
    plt.xlabel('real(lambda)')
    plt.ylabel('imag(lambda)')
    
    filename='plot({0}=N)({1}=p)({2}=k)({3}=mat).png'.format(N,p,k,mat)
    plt.savefig(filename, bbox_inches='tight')
    plt.close('all')
    # print("plotted")
    return

def solveD(gamma,p,N,k,b,runs,skip=False):
    reD=[]
    imD=[]
    for i in range(runs):
        if (i%10==0):
            print i/1000.0 
        eigvalD=np.linalg.eigvals(D(gamma,p,N,k))
        for z in eigvalD:
            if(skip and z.real>3):
                print('skipped in D')
            else:
                reD.append(z.real)
                imD.append(z.imag)
    plot(reD,imD,b,N,p,k,'D')

def solveA(gamma,p,N,k,b,runs,skip=False):
    reA=[]
    imA=[]
    for i in range(runs):
        if (i%10==0):
            print i/1000.0 
        eigvalA=np.linalg.eigvals(A(p,N))
        for z in eigvalA:
            if(skip and z.real>3):
                print('skipped in A')
            else:
                reA.append(z.real)
                imA.append(z.imag)
    plot(reA,imA,b,N,p,k,'A')
    
def solveS(gamma,p,N,k,b,runs,skip=False):
    reS=[]
    imS=[]
    for i in range(runs):
        if (i%10==0):
            print i/1000.0 
        eigvalS=np.linalg.eigvals(S(gamma,p,N,k))
        for z in eigvalS:
            if(skip and z.real>3):
                print('skipped in S')
            else:
                reS.append(z.real)
                imS.append(z.imag)
    plot(reS,imS,b,N,p,k,'S')
    
def solveDequ(p,N,k,b,runs,skip=False): 
    reDequ=[]
    imDequ=[]  
    x=0
    if skip:
        x=1
    for i in range(x,N): #skipping number 0
        z=lmbdaD(i,k,N,p)
        reDequ.append(z.real)
        imDequ.append(z.imag)
    plot(reDequ,imDequ,b,N,p,k,'Dequ')


"""
k is connectivity its even
0<p<1 
N is size of the matricies
b is bins
"""  
k=4
p=.2 
N=1000
gamma=1
b=20
runs=1

#for j in range(2,1000,100):
 #   for i in range(4):
solveD(gamma,p,N,k,b,runs)
solveA(gamma,p,N,k,b,runs)
solveS(gamma,p,N,k,b,runs)
solveDequ(p,N,k,b,runs)


