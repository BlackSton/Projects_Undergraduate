from __future__ import print_function,division
from random import random
from mpl_toolkits.mplot3d import Axes3D
from math import sin,cos,pi,acos
import matplotlib.pyplot as plt
from numpy import arange
from pylab import plot,show

N=1500
dA=0.01
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x=[]
y=[]
z=[]
for i in range(N):
    phi = 2*pi*random()                 # acos(1-2*random())
    psi = acos(1-2*random())            # °¡ÁßÄ¡
    x.append(sin(psi)*cos(phi))
    y.append(sin(psi)*sin(phi))
    z.append(cos(psi))

X=[]
Y=[]
A=-1.0
a = 0
y_sum = 0
y_avg = []
for A in arange(-1,1,dA):
    a +=1
    count=0
    X.append(A)
    for i in range(N):
        if z[i]>=A and z[i]<A+dA:
            count+=1
    y_sum += count
    Y.append(count)
    y_avg.append(y_sum/a)

ax.scatter(x,y,z, c='r', marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plot(X,Y,X,y_avg)
show()
