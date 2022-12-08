from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random,randint
from math import acos,sqrt

def distance(pos1,pos2):
    return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)

def V_LJ_a(pos1,pos2):
    global epsilon,r_c,sigma,Matom
    a = array([0.,0.,0.])
    r_d = distance(pos1,pos2)
    r = pos1 - pos2
    a = ((24*epsilon*r/sigma**2)*(2*(sigma/r_d)**14-(sigma/r_d)**8)\
        -(24*epsilon*r/sigma**2)*(2*(sigma/r_c)**14-(sigma/r_c)**8))/Matom
    return a
def Verlet(pos1,pos2,v1,v2,a):
    global dt,r_c,sigma
    if distance(pos1,pos2) < r_c and sigma*0.8 < distance(pos1,pos2):
        dx1 = (v1 * dt) + (0.5 * a(pos1,pos2) * dt**2)
        dx2 = (v2 * dt) + (0.5 * a(pos2,pos1) * dt**2)
        V_a_dx1 = a(pos1+dx1,pos2+dx2)
        V_a_dx2 = a(pos2+dx2,pos1+dx1)
        dv1 = 0.5*(a(pos1,pos2) + V_a_dx1)*dt
        dv2 = 0.5*(a(pos2,pos1) + V_a_dx2)*dt
    else:
        dx1 = v1 * dt
        dx2 = v2 * dt
        dv1,dv2 = 0.,0.
    return dx1,dx2,dv1,dv2

winsize =500 #â ������ ����
epsilon = 10.2
sigma = 0.6
r_c = 2.5*sigma
# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
natom = 4 # �� ���� ���ڰ���
Natom = natom**3 #�ѿ��� ����
Ratom = 0.6  #���� ������(����:m)
Matom = 1. # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
k = 1. # ������ ���
dt = 1.E-4 #�ð� ����
T =100. #�µ�(300K)
Lc=1.2
## �ùķ��̼� â���� ##
scene = display(title="�з¿� ���� ����", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))
deltav = 1.
vdist = gdisplay(x=0, y=winsize, ymax = Natom*deltav/10.,
             width=winsize, height=0.6*winsize, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan)

dv = 1.
for v in arange(0.,50.+dv,dv): # theoretical prediction
    theory.plot(pos=(v,
        (deltav/dv)*Natom*4.*pi*((Matom/(2.*pi*k*T))**1.5)
                     *exp((-0.5*Matom*v**2)/(k*T))*v**2*dv))

observation = ghistogram(bins=arange(0.,50.,deltav),
                        accumulate=1, average=1, color=color.red)

## �ܺ� �� ���� ##
xaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L-L/Lc,L/Lc,L/Lc)], color=color.blue)
yaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L/Lc,L-L/Lc,L/Lc)], color=color.blue)
zaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L/Lc,L/Lc,L-L/Lc)], color=color.blue)
xaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L/Lc,L-L/Lc,L-L/Lc), (L/Lc,L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L-L/Lc)], color=color.blue)
yaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L/Lc), (L-L/Lc,L-L/Lc,L/Lc)], color=color.blue)
zaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L-L/Lc,L-L/Lc,L/Lc), (L/Lc,L-L/Lc,L/Lc), (L/Lc,L-L/Lc,L-L/Lc)], color=color.blue)
    
atom = [] #���� 
atom_Pos = [] #���� ��ġ
atom_p = [] # ���� ���

## ������ ��ġ�� �ӵ� ����##
## ������ ��ġ ���� ���� ����##
for i in range(natom):
    for j in range(natom):
        for K in range(natom):
            x = ((L-7*Ratom)/(natom-1))*i+3.5*Ratom
            y = ((L-7*Ratom)/(natom-1))*j+3.5*Ratom
            z = ((L-7*Ratom)/(natom-1))*K+3.5*Ratom
            atom_Pos.append((x,y,z))
            atom.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))
            ## �ӵ��� ������ ���� ����ȭ ##
            phi = random()*2*pi
            psi = acos(1-2*random())        # ������
            pavg = sqrt(2.*Matom*1.5*k*T) #������ ��� �ӵ� ����
            atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

Pos = array(atom_Pos)
v = array(atom_p)/Matom
while True:
    observation.plot(data=mag(v))
    rate (1000)


    ## �� ��� �浹 ��� ##
    outside = greater_equal(Pos,L-Ratom)
    v1 = v*outside
    v = v-v1-abs(v1)
    ## �عٴ� ��� �浹 ��� ##
    outside = less_equal(Pos,Ratom)
    v1 = v*outside
    v = v-v1+abs(v1)
    ## ���� ��ġ �缳�� ##   
    for i in range(Natom):
        atom[i].pos = Pos[i]
