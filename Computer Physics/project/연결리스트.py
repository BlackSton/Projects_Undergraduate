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
sigma = 1
r_c = 2.5*sigma
# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
natom = 5 # �� ���� ���ڰ���
Natom = natom**3 #�ѿ��� ����
Ratom = 0.6  #���� ������(����:m)
Matom = 1. # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
k = 1. # ������ ���
dt = 1.E-4 #�ð� ����
T =100. #�µ�(300K)
## �ùķ��̼� â���� ##


atom = [] #���� 
atom_Pos = [] #���� ��ġ
atom_p = [] # ���� ���

## ������ ��ġ�� �ӵ� ����##
## ������ ��ġ ���� ���� ����##
for i in range(natom):
    for j in range(natom):
        for K in range(natom):
            x = ((L-2*Ratom)/(natom-1))*i+Ratom
            y = ((L-2*Ratom)/(natom-1))*j+Ratom
            z = ((L-2*Ratom)/(natom-1))*K+Ratom
            atom_Pos.append((x,y,z))

            ## �ӵ��� ������ ���� ����ȭ ##
            phi = random()*2*pi
            psi = acos(1-2*random())        # ������
            pavg = sqrt(2.*Matom*1.5*k*T) #������ ��� �ӵ� ����
            atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

Pos = array(atom_Pos)       # �� ���� ��ġ ���
v = array(atom_p)/Matom     # �� ���� �ӵ� ���


Lc=int(L//r_c)                   # ĭ�� ����    
head=zeros(Lc**3)      # head ĭ�����
mc=zeros(3,int)                 # ĭ ��ġ
lscl=zeros(Natom)           # lscl ĭ�����

for c in range(Lc**3):
    head[c]=-1
for i in range(Natom):
    for a in range(3):
        mc[a]=floor(Pos[i,a]/(L/Lc))      # ��� �κп��� ���� �����ϰ� �ϱ����عٲ�
    c=mc[0]*Lc**2+mc[1]*Lc+mc[2]
    # print(c,Pos[i,a],mc[0],mc[1],mc[2],Pos[i,0],Pos[i,0]/r_c)          # �� ���ڴ� ��� �������� ��Ÿ��
    lscl[i]=head[c]                       # �ڿ� head �޾ƸԱ�
    head[c]=i                             # �� �밡�� ���ö����� �ֽ�ȭ
    
print(head)

print(lscl)

while True:
    for i in range(Natom):
        for j in arange(i+1,Natom,1):
            dx0,dx1,dv0,dv1 = Verlet(Pos[i],Pos[j],v[i],v[j],V_LJ_a)
            v[i] = v[i] + dv0
            v[j] = v[j] + dv1
            Pos[i] = Pos[i] + dx0
            Pos[j] = Pos[j] + dx1
