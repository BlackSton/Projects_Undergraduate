from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

winsize =500 #â ������ ����

# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����(����:m)
B_L = 0.5 # ���� �α�(����:m)
Natom = 30 #���� ����
Ratom = 0.1  #���� ������(����:m)
Matom = 4E-3/6E23 # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
g = 9.8 # �߷»��
k = 1.4E-23 # ������ ���
dt = 1.E-4 #�ð� ����
beta =1. # �������
T =300 #�µ�(300K)

# �׷��� ���� ���� ����
time=0 # �ʱ� �ð�
theory_sum = 0 #��� �з�������
theory_avg = [] #��� �з�-�׷���-
M_B_Min = 1E-20 #�ּ� ����
M_B_Max = 5E-20 #�ְ� ����
dM = 5E-21 #���� ����
graph_color =[] #�׷��� ����

## �׷��� ���� �� ##
Time_to_P = gdisplay(title="�ð��� ���� �з�",x = 0,y = winsize,ymax = M_B_Max/500,\
                  width = winsize, height = 0.6*winsize, xtitle ="Time",ytitle="P") #�ð��� ���� �з� �׷��� ȭ�� ����
Time_P = gdots(color=color.red)#�ð��� ���� �з� �׷��� �� ��°�
T_to_P = gdisplay(title="�з°� ���ǰ��� ����",x = winsize,y = 0,ymax = L**3,\
                  width = winsize, height = winsize, xtitle ="P",ytitle="V") #�µ��� �з°���  �׷���ȭ�� ����
P_T_Theory = gcurve(color = color.cyan) #�̷� �׷��� ���� ����
P_T_Experimental = gdots(color=color.red) #���� �׷��� ���� ����
for graph_M in arange(dM,M_B_Max+1E-22+dM,dM): # �̷� ���� �׷��� ����
    P_T_Theory.plot(pos=(graph_M,Natom*k*L**3/(graph_M*9.8)))

time_to_P = gdisplay(title="�ð��� ���� ����",x=winsize, y=winsize,\
                     ymax = L**3,width=winsize, height=0.6*winsize, xtitle='Time', ytitle='V') #�ð��� ���� ���� �׷��� ȭ�� ����
theory = gdots(color=color.red)
for i in range(int((M_B_Max-M_B_Min+dM)/dM)):
    graph_color.append((0.5,i/8.,i/8.+0.3))
    theory_avg.append(gcurve(color =graph_color[i]))

for M in arange(M_B_Min,M_B_Max+1E-22,dM):  #M_B_Min ���� B_M_Max ���� dM �������� �ùķ��̼�..
    ## �ùķ��̼� â���� ##
    scene = display(title="�з¿� ���� ����", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

    ## �ܺ� �� ���� ##
    xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.blue)
    yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.blue)
    zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.blue)
    xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.blue)
    yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.blue)
    zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.blue)

    ## õ�� ���� ##
    T_B = box(pos=vector(L/2,L,L/2),size=vector(L,B_L,L))
    B_Pos = vector(L/2,L,L/2)
    B_v = vector(0,0,0)
    
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
    p = array(atom_p)
    theory_length = 500 # ����� ���� ����
    theory_P_avg = empty([theory_length],float)# ��������� �� ����
    theory_P_avg[:]=1E-23# ��� �ʱⰪ(�µ������� ��������� ��.<-ù ��վз°��� �׷��� �������� �̻��ϸ� �̰� �����ϸ� ��.)
    while time < 3.0:
        time+= dt
        rate (1000)
        ## dt�� �ð��� ������ ��ġ ��� ##
        Pos = Pos + (p/Matom)*dt
        B_v = B_v - vector(0,(9.8)*dt,0)
        ## �� ��� �浹 ��� ##
        outside = greater_equal(Pos,[L-Ratom,B_Pos[1]-Ratom-B_L,L-Ratom])
        p1 = p*outside
        p = p-p1-abs(p1)
        B_V = B_v[1] #�浹�� ���� �ӵ�
        theory_P = 0
        for c in range(Natom):
            if outside[c][1]:
                b_v = -p[c][1]/Matom # �浹�� ���� �ӵ�
                theory_P += b_v*Matom
                b_v0= (b_v*(Matom-M)+2*M*B_V)/(Matom+M)#�浹�� ���� �ӵ�
                B_V0= (B_V*(M-Matom)+2*Matom*b_v)/(Matom+M)#�浹�� ���� �ӵ�
                B_V = B_V0
                theory_P += b_v *Matom
                p[c][1] = b_v0*Matom
        B_v = vector(0,B_V,0)
        B_Pos = B_Pos + B_v*dt
        ## �ð��� ���� �з� �� ���(������ �з¸� ����ϴ°���)##
        theory_P_avg=append(theory_P_avg,theory_P)
        theory_P_avg=delete(theory_P_avg,0)
        Time_P.plot(pos=(time,sum(theory_P_avg)/theory_length),color=graph_color[int((M-M_B_Min)/dM)])
        ## �عٴ� ��� �浹 ��� ##
        outside = less_equal(Pos,Ratom)
        p1 = p*outside
        p = p-p1+abs(p1)
        ## �ð��� ���� ���� �� ��� ##
        theory_volume= L*L*(B_Pos[1]-B_L)
        theory_avg[int((M-M_B_Min)/dM)].plot(pos=(time,theory_volume))
        ## ���� ��ġ �缳�� ##   
        for i in range(Natom):
            atom[i].pos = Pos[i]
        T_B.pos = B_Pos
    P_T_Experimental.plot(pos=(M,theory_volume),color=graph_color[int((M-M_B_Min)/dM)]) #������ �׷���ȭ

    ## �� �ʱ�ȭ ##
    time = 0
    theory_sum = 0
    scene.delete()

while 1:
    rate(1)
