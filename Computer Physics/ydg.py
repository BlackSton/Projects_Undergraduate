from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random
from math import acos,sqrt

def distance(pos1,pos2):
    return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)
def c1_select(c): # c1 �� ��ȯ ���ִ� �Լ�
    global Lc
    a_s = [] # c1�� �迭
    for i in range(3):
        a_z = -Lc**2 + i*Lc**2 # Z�� +- ��Ű�°�
        a = arange(a_z+c-Lc-1,a_z+c-Lc+2) # x�ప ����
        a1 = append(a,a+Lc) # y�� +
        a_2 = append(a1,a+Lc*2) #y�� ++
        a_s = append(a_s,a_2) #z�࿡ �߰�
    return a_s
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
sigma = 0.5
r_c = 2.5*sigma
# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
natom = 5 # �� ���� ���ڰ���
Natom = natom**3 #�ѿ��� ����
Ratom = 0.5  #���� ������(����:m)
Matom = 1. # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
k = 1. # ������ ���
dt = 1.E-4 #�ð� ����
T =100. #�µ�(300K)
## �ùķ��̼� â���� ##
scene = display(title="�з¿� ���� ����", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.blue)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.blue)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.blue)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.blue)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.blue)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.blue)

atom = [] #���� 
atom_Pos = [] #���� ��ġ
atom_p = [] # ���� ���

## ������ ��ġ�� �ӵ� ����##
## ������ ��ġ ���� ���� ����##
for i in range(natom):
    for j in range(natom):
        for K in range(natom):
            x = ((L-7*Ratom)/(natom-1))*i+3.5*Ratom # ������ ������ 3.5*���ڹ��������� �����.
            y = ((L-7*Ratom)/(natom-1))*j+3.5*Ratom # ������ ������ 3.5*���ڹ��������� �����.
            z = ((L-7*Ratom)/(natom-1))*K+3.5*Ratom # ������ ������ 3.5*���ڹ��������� �����.
            atom_Pos.append((x,y,z))
            atom.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))
            ## �ӵ��� ������ ���� ����ȭ ##
            phi = random()*2*pi
            psi = acos(1-2*random())        # ������
            pavg = sqrt(2.*Matom*1.5*k*T) #������ ��� �ӵ� ����
            atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

Pos = array(atom_Pos)       # �� ���� ��ġ ���
v = array(atom_p)/Matom     # �� ���� �ӵ� ���


Lc=int(2+L//r_c)            # ���࿡���� ĭ�� ����    +2 �ϴ� ������ �� �ܺ����ο����� ������ֱ� ���ؼ�.
mc=empty(3,int)             # ���࿡���� ĭ ��ġ
lscl=zeros(Natom,int)           # lscl ĭ�����
head=empty(Lc**3,int) 

while True:
    rate(1000)
    #���Ḯ��Ʈ ����#
    for c in range(Lc**3):
        head[c]=-1
    for i in range(Natom):
        for a in range(3):
            mc[a]=floor(Pos[i,a]/(L/Lc))      # ��� �κп��� ���� �����ϰ� �ϱ����عٲ�
        c=mc[0]*Lc**2+mc[1]*Lc+mc[2]
        lscl[i]=head[c]                       # �ڿ� head �޾ƸԱ�
        head[c]=i
        
    # ��ȣ �ۿ� ���� �ʴ� ���ڵ��� �׳� �����δ�.
    for cc in range(Natom):
        if lscl[cc] == -1:
            Pos[cc] = Pos[cc] + v[cc]*dt
            
    #��ȣ�ۿ� ��� �ڵ�#
    for c in range(Lc**3):
        i = head[c]
        while i!= -1:
            for c1 in c1_select(c):
                j = head[int(c1)]
                while j != -1:
                    if i < j:
                        if distance(Pos[i],Pos[j])<r_c:
                            #��ȣ�ۿ�
                            dx0,dx1,dv0,dv1 = Verlet(Pos[i],Pos[j],v[i],v[j],V_LJ_a)
                            v[i] = v[i] + dv0
                            v[j] = v[j] + dv1
                            Pos[i] = Pos[i] + dx0
                            Pos[j] = Pos[j] + dx1
                    j = lscl[j]
            i = lscl[i]
    
        ## �� ��� �浹 ��� ##
    outside = greater_equal(Pos,L-Ratom-L/Lc) # L/Lc�� ��ĭ�� ����
    v1 = v*outside
    v = v-v1-abs(v1)
    ## �عٴ� ��� �浹 ��� ##
    outside = less_equal(Pos,Ratom+L/Lc) # L/Lc�� ��ĭ�� ����
    v1 = v*outside
    v = v-v1+abs(v1)
    for i in range(Natom):
        atom[i].pos = Pos[i]
