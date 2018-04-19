from math import *
from numpy import *

N=2
#w1=1.12539373e+12
#w2=9.51411015e+11
#zet1=2.68373170e-05
#zet2=5.62391553e-05


m=zeros(2*N)
#m[0]=w1+w2
#m[1]=(w1*zet1+w2*zet2)/m[0]
#m[2]=(w1*zet1**2+w2*zet2**2)/m[0]
#m[3]=(w1*zet1**3+w2*zet2**3)/m[0]

#m[0]=2.0768047422e+12
#m[1]=0.1/m[0]
#m[2]=9.630178318e-15/m[0]
#m[3]=1.391105016e-27/m[0]
#m[4]=2.679317874062885e-40/m[0]
#m[5]=6.450577224762668e-53/m[0]
#m[6]=1.863606267942968e-65/m[0]
#m[7]=6.281401236488749e-78/m[0]

#m[0] = 25.8569547
#m[1] = 0.0489969/m[0]
#m[2] = 9.549288328E-5/m[0]
#m[3] = 1.9098593171E-7/m[0]

m[0] = 1.0
m[1] = 1.0
m[2] = 1.0
m[3] = 1.0

print m

Moments=zeros(2*N)
Moments[0] = 1.0
for i in xrange(1,2*N):
    Moments[i] = (-1)**i * m[i]

B=zeros((2*N, 2*N+1))
B[0,0]=1.0
B[:,1] = Moments[:]
for j in range(2,2*N+1):
    for i in range(2*N+1-j):
        B[i,j] = B[0,j-1]*B[i+1,j-2] - B[0,j-2]*B[i+1,j-1]

c=zeros(2*N)
c[0]=0.0
for i in xrange(1,2*N):
    c[i] = B[0,i+1]/(B[0,i]*B[0,i-1])

A=zeros((N,N))
for i in range(0,N-1):
    A[i,i] = c[2*i+1] + c[2*i]
    A[i,i+1] = sqrt(abs(c[2*i+2]*c[2*i+1]))
    A[i+1,i] = sqrt(abs(c[2*i+2]*c[2*i+1]))
i=N-1
A[i,i] = c[2*i+1] + c[2*i]

for i in range(N):
   print '\n'
   for j in range(N):
      print A[i,j], ' ',

evalues, evectors = linalg.eig(A)
print evalues 
print evectors
weight=zeros(N)
abscissa=zeros(N)
abscissa=evalues
for i in range(N):
    weight[i]=m[0]*evectors[0,i]**2
print '\n'
print weight,abscissa,weight*abscissa

# Checking values of matirx A on inversion
import numpy as np 
x1 = abscissa[0]
x2 = abscissa[1]
A = np.matrix([[1,1,0,0],[0,0,1,1],[-(x1)**2,-(x2)**2,2*x1,2*x2],[-(x1)**3,-(x2)**3,3*x1,3*x2]])
inverse_A = np.linalg.inv(A)
print inverse_A
