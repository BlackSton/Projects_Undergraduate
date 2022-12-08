from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

winsize =500 #â ������ ����

# �Ϲ� ��� ����
L = 10 # ��� �Ѻ� ����
Natom = 30 #���� ����
Ratom = 0.5  #���� ������
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
deltav = 100. # �ӵ� ���� �׷��� �׸���
vdist = gdisplay(x=0, y=winsize, ymax = Natom*deltav/1000.,width=winsize, height=0.6*winsize, xtitle='v', ytitle='dN')
dN_v = gcurve(color=color.cyan)
dv = 10.
observation = ghistogram(bins=arange(0.,3000.,deltav),accumulate=1, average=1, color=color.red)
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
    for v in arange(0.,3001.+dv,dv): # �̷� ���� �׷���
        dN_v.plot(pos=(v,(deltav/dv)*Natom*4.*pi*((Matom/(2.*pi*k*T))**1.5)*exp((-0.5*Matom*v**2)/(k*T))*v**2*dv))

    ## �ܺ� �� ���� ##
    xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.orange)
    yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.orange)
    zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.orange)
    xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.orange)
    yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.orange)
    zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.orange)

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

    while time < 0.2:
        time+= dt
        rate (100)
        observation.plot(data=mag(p/Matom))
        ## dt�� �ð��� ������ ��ġ ��� ##
        Pos = Pos + (p/Matom)*dt

        ## ���ڵ��� �浹 ���� �Ǵ� ##
        atom_to_atom_Pos = Pos-Pos[:,newaxis] # �����ڰ��� ��ǥ(x,y,z) �Ÿ�
        atom_to_atom_distance = sqrt(sum(square(atom_to_atom_Pos),-1)) # �������� ���� �Ÿ�
        hit = less_equal(atom_to_atom_distance,2.*Ratom)-identity(Natom) # ���ڵ鰣�� �浹���� ����(�ο��ڻ����� ���� < �������� ������ ������ ��)
        hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*Natoms+j
        ## ���ڵ��� �浹 ��� ##
        for ij in hitlist:
            i, j = divmod(ij,Natom) # ������ �з�
            hitlist.remove(j*Natom+i) # �浹����Ʈ���� ����
            ptot = p[i]+p[j] #�浹�ϴ� �� ������ ��� ��
            #������ �ӵ� ���#
            vi = p[i]/Matom
            vj = p[j]/Matom
            #���ڵ��� �ε����� ������ ��ġ�� ��ȯ#
            a = mag(vj-vi)**2
            b = 2*dot(Pos[i]-Pos[j],vj-vi)
            c = mag(Pos[i]-Pos[j])**2-(2.*Ratom)**2
            d = b**2-4.*a*c
            if d < 0: continue # something wrong; ignore this rare case
            deltat = (-b+sqrt(d))/(2.*a)
            Pos[i] = Pos[i]-(p[i]/Matom)*deltat
            Pos[j] = Pos[j]-(p[j]/Matom)*deltat
            mtot = Matom * 2. # �� ������ �� ����
            #���� �߽� �������� ��� ���#
            pcmi = p[i]-ptot*Matom/mtot
            pcmj = p[j]-ptot*Matom/mtot
            rrel = norm(Pos[j]-Pos[i]) #���ڻ��ǰŸ� �������ͷ� ����
            #�浹�� �����߽ɱ��� ��� ��ȯ
            pcmi = pcmi-2*dot(pcmi,rrel)*rrel 
            pcmj = pcmj-2*dot(pcmj,rrel)*rrel
            #������ ������ ��� ��ȯ
            p[i] = pcmi+ptot*Matom/mtot 
            p[j] = pcmj+ptot*Matom/mtot
            #�浿�� ��ġ ����
            Pos[i] = Pos[i]+(p[i]/Matom)*deltat
            Pos[j] = Pos[j]+(p[j]/Matom)*deltat

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
