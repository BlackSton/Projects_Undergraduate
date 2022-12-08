from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

winsize =500 #â ������ ����

# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����
Natom = 30 #���� ����
Ratom = 0.4  #���� ����
Matom = 4E-3/6E23 # ���� ���� 4g/�ƺ�����μ� = 0.004/6.02*10^23
k = 1.4E-23 # ������ ���
dt = 1.E-4 #�ð� ����

# �׷��� ���� ���� ����
time=0 # �ʱ� �ð�
theory_sum = 0 #��� �з�������
theory_avg = [] #��� �з�-�׷���-
T_Min = 50 #�ּ� �µ�
T_Max = 500 #�ְ� �µ�
dT = 50. #�µ� ����
graph_color =[] #�׷��� ����

## �׷��� ���� �� ##
T_to_P = gdisplay(title="�µ��� �з°��� ����",x = winsize,y = 0,ymax = T_Max*Natom*k/L**3,\
                  width = winsize, height = winsize, xtitle ="T",ytitle="P") #�µ��� �з°���  �׷���ȭ�� ����
P_T_Theory = gcurve(color = color.cyan) #�̷� �׷��� ���� ����
P_T_Experimental = gdots(color=color.red) #���� �׷��� ���� ����
for graph_T in arange(0.,501.+dT,dT): # �̷� ���� �׷��� ����
    P_T_Theory.plot(pos=(graph_T,graph_T*Natom*k/L**3))

time_to_P = gdisplay(title="�ð��� ���� ��� �з�",x=winsize, y=winsize,\
                     ymax = T_Max*Natom*k/L**3,width=winsize, height=0.6*winsize, xtitle='Time', ytitle='p') #�ð��� ���� ��� �з� �׷���ȭ�� ����
theory = gdots(color=color.red)
for i in range(int((T_Max-T_Min+dT)/dT)):
    graph_color.append((0.5,i/8.,i/8.+0.3))
    theory_avg.append(gcurve(color =graph_color[i]))

for T in arange(T_Min,T_Max+1,dT):  #T_Min ���� T_Max ���� dT �������� �ùķ��̼�..

    ## �ùķ��̼� â���� ##
    scene = display(title="�µ��� ���� �з�", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

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
    p = array(atom_p)

    while time < 0.04:
        time+= dt
        rate (100)
        ## dt�� �ð��� ������ ��ġ ��� ##
        Pos = Pos + (p/Matom)*dt

        ## �عٴ� ��� �浹 ��� ##
        outside = less_equal(Pos,Ratom)
        p1 = p*outside
        p = p-p1+abs(p1)

        ## �عٴ� �浹�� �з� ���� �׷��� ��� ##
        theory_sum += (sum(abs(2*p*outside)))/dt
        theory_sum_avg = theory_sum/((time/dt)*3*L**2)
        #theory.plot(pos=(time,graph/(dt*3*L**2)))
        theory_avg[int((T-T_Min)/dT)].plot(pos=(time,theory_sum_avg))

        ## �� ��� �浹 ��� ##
        outside = greater_equal(Pos,L-Ratom)
        p1 = p*outside
        p = p-p1-abs(p1)
        
        ## ���� ��ġ �缳�� ##   
        for i in range(Natom):
            atom[i].pos = Pos[i]
   
    print("�µ�:",T,"���� ��� �з�:",theory_sum_avg,"�̷� ��� �з�:",T*Natom*k/L**3,"������",100*abs(theory_sum_avg-T*Natom*k/L**3)/(T*Natom*k/L**3),"%")

    P_T_Experimental.plot(pos=(T,theory_sum_avg),color=graph_color[int((T-T_Min)/dT)]) #������ �׷���ȭ

    ## �� �ʱ�ȭ ##
    time = 0
    theory_sum = 0
    scene.delete()

while 1:
    rate(1)
