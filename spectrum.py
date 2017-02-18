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

def D(gamma,p,N,k):
    d= (gamma*C(k,N)+alpha(p)*J(N))/np.sqrt(N)
    #print(d)
    return d
    
def S(gamma,p,N,k):
    s=A(p,N)+D(gamma,p,N,k)
   # print(s)
    return s

def plot(re,im,b,N,p,k,mat):
    plt.figure(1)
    plt.subplot(211)
    plt.hist2d(re, im, bins=b,norm=LogNorm())
    plt.colorbar()
    title='Eigenvalues of {3},N={0},p={1},k={2},{4}of them'.format(N,p,k,mat,len(re))
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

    
def solve(gamma,p,N,k,b,runs,mat,skip=False):
    re=[]
    im=[]
    for i in range(runs):
        #timer
        #if (i%10==0):
         #   print i/ float(runs)
        #choose your set of eigenvalues
        eigval=[]
        if mat=='S':
            eigval=np.linalg.eigvals(S(gamma,p,N,k))
        elif mat=='A':
            eigval=np.linalg.eigvals(A(p,N))
        elif mat=='D':
            eigval=np.linalg.eigvals(D(gamma,p,N,k))
        elif mat=='C':
            eigval=np.linalg.eigvals(C(k,N))
        elif mat=='J':
            eigval=np.linalg.eigvals(J(N))
        
        else:
            print("need a matrix type")
            
        for z in eigval:
            if(skip and z.real>3):
                print('skipped in '+mat)
            else:
                re.append(z.real)
                im.append(z.imag)
    plot(re,im,b,N,p,k,mat)
    
def solveDequ(p,N,k,b,runs,skip=False): 
    reDequ=[]
    imDequ=[]  
    x=0
    if skip:
        x=1
        print('skipped in Dequ')
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
k=40
p=.5 
N=100
gamma=1
b=20
runs=100

solve(gamma,p,N,k,b,runs,'D')
        #print("D")
solve(gamma,p,N,k,b,runs,'C')
#print("C")
solve(gamma,p,N,k,b,runs,'J')
        #print("J")
solve(gamma,p,N,k,b,runs,'A')
        #print("A")
solve(gamma,p,N,k,b,runs,'S')
        #print("S")
solveDequ(p,N,k,b,runs)
        #print("Dequ")



