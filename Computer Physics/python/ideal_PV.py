from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

def distance(pos1,pos2):
    return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)

def V_LJ_a(pos1,pos2):
    global epsilon,r_c,sigma,Matom
    a = array([0.,0.,0.])
    r_d = distance(pos1,pos2)
    if r_d < r_c:
        r = pos1 - pos2
        a = ((24*epsilon*r/sigma**2)*(2*(sigma/r_d)**14-(sigma/r_d)**8)\
             -(24*epsilon*r/sigma**2)*(2*(sigma/r_c)**14-(sigma/r_c)**8))/Matom
    else:
        a = array([0.,0.,0.])
    return a
def Verlet(pos1,pos2,v1,v2,a):
    global dt
    dx1 = (v1 * dt) + (0.5 * a(pos1,pos2) * dt**2)
    dx2 = (v2 * dt) + (0.5 * a(pos2,pos1) * dt**2)
    V_a_dx1 = a(pos1+dx1,pos2+dx2)
    V_a_dx2 = a(pos2+dx2,pos1+dx1)
    dv1 = 0.5*(a(pos1,pos2) + V_a_dx1)*dt
    dv2 = 0.5*(a(pos2,pos1) + V_a_dx2)*dt
    return dx1,dx2,dv1,dv2

winsize =500 #â ������ ����
epsilon = 1.2
sigma = 0.5
r_c = 2.5*sigma
# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
Natom = 3 #���� ����
Ratom = 0.1  #���� ������(����:m)
Matom = 4E-3/6E23 # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
g = 9.8 # �߷»��
k = 1.4E-23 # ������ ���
dt = 1.E-4 #�ð� ����
T =3 #�µ�(300K)

## �ùķ��̼� â���� ##
scene = display(title="�з¿� ���� ����", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

## �ܺ� �� ���� ##
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
for i in range(Natom):
    ## ������ ��ġ ����ȭ ##
    x = Ratom+random()*(L-2*Ratom)
    y = Ratom+random()*(L-2*Ratom)
    z = Ratom+random()*(L-2*Ratom)
    atom_Pos.append((x,y,z))
    atom.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))

    ## �ӵ��� ������ ���� ����ȭ ##
    phi = random()*2*pi
    psi = random()*2*pi
    pavg =sqrt(3.*Matom*k*T) #������ ��� �ӵ� ����
    atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

Pos = array(atom_Pos)
v = array(atom_p)/Matom
while True:
    rate (1000)
    Pos = Pos + v*dt
    r = Pos-Pos[:,newaxis] # all pairs of atom-to-atom vectors
    rmag = sqrt(sum(square(r),-1)) # atom-to-atom scalar distances
    hit = less_equal(rmag,r_c)-identity(Natom)    #radius+radius[:,None]
    hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*Natoms+j
    for ij in hitlist:
        i, j = divmod(ij,Natom) # decode atom pair
        hitlist.remove(j*Natom+i) # remove symmetric j,i pair from list
        dx0,dx1,dv0,dv1 = Verlet(Pos[i],Pos[j],v[i],v[j],V_LJ_a)
        v[i] = v[i] + dv0
        v[j] = v[j] + dv1
        Pos[i] = Pos[i] + dx0
        Pos[j] = Pos[j] + dx1        
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
while 1:
    rate(1)
