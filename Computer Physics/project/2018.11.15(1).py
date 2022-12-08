from __future__ import print_function,division
from random import random
from mpl_toolkits.mplot3d import Axes3D
from math import sin,cos,pi
import matplotlib.pyplot as plt
from numpy import arange
from pylab import plot,show

N=1500
dA=0.05
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # 3차원 그래프
x=[]
y=[]
z=[]
for i in range(N):
    phi = random()*2*pi                     # x,y
    psi = random()*2*pi                     # z
    x.append(sin(psi)*cos(phi))
    y.append(sin(psi)*sin(phi))
    z.append(cos(psi))

X=[]
Y=[]
A=-1.0
for A in arange(-1,1,dA):  #0.9정도로 두기(A+0.1이상일때)
    count=0
    X.append(A)
    for i in range(N):
        if z[i]>=A and z[i]<A+dA:
            count+=1
    Y.append(count)

ax.scatter(x,y,z, c='r', marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plot(X,Y)
show()
