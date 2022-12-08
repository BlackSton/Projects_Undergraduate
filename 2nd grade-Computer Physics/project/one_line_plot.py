from pylab import plot,show
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
Ratom =0.5 # ���� ������
Matom = 1. # ���� ����
epsilon = 10.2
sigma = 1.
l = 1.9
r_c = 2.5*sigma
Time = 0 # �ʱ� �ð�
dt = 0.001                                   # �ð� ����
x_v = np.arange(100.,500.,10.)
y_p = []
for y in x_v:
    atom = [] # ����
    atom_pos = [] # ���� ��ġ
    atom_v = [] # ���� �ӵ�
    atom_pos.append((l,0.,0.))
    atom_pos.append((0.,0.,0.))
    atom_v.append((-y,0,0))
    atom_v.append((0.,0,0))
    pos = np.array(atom_pos)
    v = np.array(atom_v)
    first_p = abs(v[0][0])+abs(v[1][0])
    
    while Time<2*l/y:
        Time += dt
        dx0,dx1,dv0,dv1 = Verlet(pos[0],pos[1],v[0],v[1],V_LJ_a)
        v[0] = v[0] + dv0
        v[1] = v[1] + dv1
        pos[0] = pos[0] + dx0
        pos[1] = pos[1] + dx1
    last_p = abs(v[0][0])+abs(v[1][0])
    y_p.append(last_p/first_p)
    Time = 0
plot(x_v,y_p)
show()
