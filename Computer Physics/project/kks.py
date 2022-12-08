from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random
from math import acos,sqrt
def counting(pos1,pos2,v1,v2):    
    s_v = 100 # �Ӱ�ӵ� <- �̼ӵ� �̻�Ѿ�� ī���ð����� �ø�.
    r_v = abs(dot(pos1-pos2,v1-v2)/sqrt(pos1[0]**2+pos1[1]**2+pos1[2]**2)) #�� ���ڿ� ���� ���ӵ� ���.
    return int(r_v/s_v)*50+1
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
def Verlet(pos1,pos2,v1,v2,a,dt):
    global r_c,sigma
    if distance(pos1,pos2) < r_c:
        dx1 = v1 * dt
        dx2 = v2 * dt
        dv1 = a(pos1+dx1,pos2+dx2)*dt
        dv2 = a(pos2+dx2,pos1+dx1)*dt
    else:
        dx1 = v1 * dt
        dx2 = v2 * dt
        dv1,dv2 = 0.,0.
    return dx1,dx2,dv1,dv2

winsize =500 #â ������ ����
epsilon = 1.2
sigma = 0.2
r_c = 2.5*sigma
# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
natom = 7 # �� ���� ���ڰ���
Natom = natom**3 #�ѿ��� ����
Ratom = 0.2  #���� ������(����:m)
Matom = 1. # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
k = 8.#8.314 # ������ ���
dt = 1.E-4 #�ð� ����
T =1000. #�µ�(300K)
## �ùķ��̼� â���� ##
scene = display(title="������ ����", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))
deltav = 5.
vdist = gdisplay(x=0, y=winsize, ymax = Natom*deltav/180.,
             width=winsize, height=0.6*winsize, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan)
dv = 1.
for v in arange(0.,800.+dv,dv): # theoretical prediction
    theory.plot(pos=(v,
        (deltav/dv)*Natom*4.*pi*((Matom/(2.*pi*k*T))**1.5)
                     *exp((-0.5*Matom*v**2)/(k*T))*v**2*dv))

observation = ghistogram(bins=arange(0.,800.,deltav),
                        accumulate=0.01, average=0.01, color=color.red)

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


Lc=int(L//r_c)            # ���࿡���� ĭ�� ����    +2 �ϴ� ������ �� �ܺ����ο����� ������ֱ� ���ؼ�.
mc=empty(3,int)             # ���࿡���� ĭ ��ġ
lscl=zeros(Natom,int)           # lscl ĭ�����
head=empty(Lc**3,int)

#�ܺ��� ����.#
xaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L-L/Lc,L/Lc,L/Lc)], color=color.blue)
yaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L/Lc,L-L/Lc,L/Lc)], color=color.blue)
zaxis = curve(pos=[(L/Lc,L/Lc,L/Lc), (L/Lc,L/Lc,L-L/Lc)], color=color.blue)
xaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L/Lc,L-L/Lc,L-L/Lc), (L/Lc,L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L-L/Lc)], color=color.blue)
yaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L-L/Lc), (L-L/Lc,L/Lc,L/Lc), (L-L/Lc,L-L/Lc,L/Lc)], color=color.blue)
zaxis2 = curve(pos=[(L-L/Lc,L-L/Lc,L-L/Lc), (L-L/Lc,L-L/Lc,L/Lc), (L/Lc,L-L/Lc,L/Lc), (L/Lc,L-L/Lc,L-L/Lc)], color=color.blue)

while True:
    observation.plot(data=mag(v))
    rate(1000)
    
    ## �� ��� �浹 ��� ##
    outside = greater_equal(Pos,L-Ratom-L/Lc) # L/Lc�� ��ĭ�� ����
    v1 = v*outside
    v = v-v1-abs(v1)
    Pos = Pos+v*dt*outside
    ## �عٴ� ��� �浹 ��� ##
    outside = less_equal(Pos,Ratom+L/Lc) # L/Lc�� ��ĭ�� ����
    v1 = v*outside
    v = v-v1+abs(v1)
    Pos = Pos+v*dt*outside
    
    #���Ḯ��Ʈ ����#
    for c in range(Lc**3):
        head[c]=-1
    for i in range(Natom):
        for a in range(3):
            mc[a]=floor(Pos[i,a]/(L/Lc))      # ��� �κп��� ���� �����ϰ� �ϱ����عٲ�
        c=mc[0]*Lc**2+mc[1]*Lc+mc[2]
        lscl[i]=head[c]                       # �ڿ� head �޾ƸԱ�
        head[c]=i
            
    #��ȣ�ۿ� ��� �ڵ�#
    for c in range(Lc**3):
        i = head[c]
        while i!= -1:
            for c1 in c1_select(c):
                j = head[int(c1)]
                while j != -1:
                    if i < j:
                        #��ȣ�ۿ�
                        if c == c1:
                            atom[i].color = color.cyan
                            atom[j].color = color.cyan
                        count = counting(Pos[i],Pos[j],v[i],v[j])
                        for kk in range(count): # verlet ����� �̿��Ͽ� ���ڰ� ��ȣ�ۿ� ���
                            dx0,dx1,dv0,dv1 = Verlet(Pos[i],Pos[j],v[i],v[j],V_LJ_a,dt/count)
                            v[i] = v[i] + dv0
                            v[j] = v[j] + dv1
                            Pos[i] = Pos[i] + dx0
                            Pos[j] = Pos[j] + dx1
                    if i == j and lscl[head[c]] == -1: # ����ȥ��������� ���.
                        Pos[i] = Pos[i] + v[i]*dt
                        atom[i].color = color.red
                    j = lscl[j]
            i = lscl[i]
    for i in range(Natom):
        atom[i].pos = Pos[i]
