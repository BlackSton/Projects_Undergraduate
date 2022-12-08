from __future__ import print_function,division
from visual import *
from visual.graph import *
from math import sqrt
import numpy as np

def V_LJ_a(pos1,pos2): # �� ���� �� ���ӵ� ���� �Լ� Lennard-Jones potential �̿�.
    global epsilon,r_c,sigma,Matom # �⺻ ��� �� ������
    a = np.array([0.,0.,0.]) # ���ӵ� ����
    r_d = sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2) # �� ���ڻ����� �Ÿ�
    if r_d < r_c: # �� ���ڰ��� �ۿ��ϴ� ���� �����̶� ������. �̶� r_c = 2.5 * sigma
        r = pos1 - pos2 # �� ��ǥ ������ ����
        a = ((24*epsilon*r/sigma**2)*((2*(sigma/r_d)**14-(sigma/r_d)**8)\
            -(2*(sigma/r_c)**14-(sigma/r_c)**8)))/Matom  # ������ �������� ���ӵ�
    else:
        a = np.array([0.,0.,0.]) # �Ÿ��� ����� �ۿ��ϴ� ���� ����.
    return a
def Verlet(pos1,pos2,v1,v2,a):
    global dt # �⺻ �� ������
    dx1 = (v1 * dt) + (0.5 * a(pos1,pos2) * dt**2) # �� ������ ��ġ ��ȭ��
    dx2 = (v2 * dt) + (0.5 * a(pos2,pos1) * dt**2) # dx = v*dt + 0.5*a(x)*dt**2
    dv1 = 0.5*(a(pos1,pos2) + a(pos1+dx1,pos2+dx2))*dt # �� ������ �ӵ� ��ȭ��
    dv2 = 0.5*(a(pos2,pos1) + a(pos2+dx2,pos1+dx1))*dt # dv = 0.5*(a(x)+a(x+dx))*dt
    return dx1,dx2,dv1,dv2 # ��ȭ�� ���
    
winsize =500 #â ������ ����
L = 10
Ratom =0.5 # ���� ������
Matom = 1. # ���� ����
epsilon = 10.2
sigma = 1.
r_c = 2.5*sigma
Time = 0 # �ʱ� �ð�
dt = 0.01 # �ð� ����
l = gdisplay(title="",x = 0,y = winsize,width = 100, height = 100) #â ������� �ʰ��ϴ¿뵵 �������
P_T_Theory = gcurve(color = color.cyan)   #â ������� �ʰ��ϴ¿뵵 �������               
for y in np.arange(0.,1.,0.1):
    scene = display(title="12-6 potential", width=winsize, height=winsize,\
                x=0, y=0,center=(0.,0.,L/2.))
    atom = [] # ����
    atom_pos = [] # ���� ��ġ
    atom_v = [] # ���� �ӵ�
    atom_pos.append((5.,y,0))
    atom_pos.append((-5.,-y,0))
    atom.append(sphere(pos=atom_pos[0],radius=Ratom,color = color.cyan\
                       ,make_trail=True,retain=1500))
    atom.append(sphere(pos=atom_pos[1],radius=Ratom,color = color.red\
                       ,make_trail=True,retain=1500))
    atom_v.append((-3.,0,0))
    atom_v.append((0.,0,0))
    pos = array(atom_pos)
    v = array(atom_v)
    print("�ʱ� ���",sqrt(v[0][0]**2+v[0][1]**2),sqrt(v[1][0]**2+v[1][1]**2))
    while Time<15/3:
        Time += dt
        rate(1000)
        dx0,dx1,dv0,dv1 = Verlet(pos[0],pos[1],v[0],v[1],V_LJ_a)
        v[0] = v[0] + dv0
        v[1] = v[1] + dv1
        pos[0] = pos[0] + dx0
        pos[1] = pos[1] + dx1
        atom[0].pos = array(pos[0])
        atom[1].pos = array(pos[1])
    print("���� ���",sqrt(v[0][0]**2+v[0][1]**2),sqrt(v[1][0]**2+v[1][1]**2))
    Time = 0
    scene.delete()
exit()
