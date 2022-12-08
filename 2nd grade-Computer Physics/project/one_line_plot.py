from pylab import plot,show
from math import sqrt
import numpy as np
def V_LJ_a(pos1,pos2): # 두 원자 의 가속도 정의 함수 Lennard-Jones potential 이용.
    global epsilon,r_c,sigma,Matom # 기본 상수 값 가져옴
    a = np.array([0.,0.,0.]) # 가속도 정의
    r_d = sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2) # 두 원자사이의 거리
    if r_d < r_c: # 두 원자간의 작용하는 힘이 조금이라도 있을때. 이때 r_c = 2.5 * sigma
        r = pos1 - pos2 # 두 좌표 사이의 벡터
        a = ((24*epsilon*r/sigma**2)*((2*(sigma/r_d)**14-(sigma/r_d)**8)\
            -(2*(sigma/r_c)**14-(sigma/r_c)**8)))/Matom  # 힘에서 질량빼면 가속도
    else:
        a = np.array([0.,0.,0.]) # 거리를 벗어나면 작용하는 힘이 없다.
    return a
def Verlet(pos1,pos2,v1,v2,a):
    global dt # 기본 값 가져옴
    dx1 = (v1 * dt) + (0.5 * a(pos1,pos2) * dt**2) # 각 원자의 위치 변화량
    dx2 = (v2 * dt) + (0.5 * a(pos2,pos1) * dt**2) # dx = v*dt + 0.5*a(x)*dt**2
    dv1 = 0.5*(a(pos1,pos2) + a(pos1+dx1,pos2+dx2))*dt # 각 원자의 속도 변화량
    dv2 = 0.5*(a(pos2,pos1) + a(pos2+dx2,pos1+dx1))*dt # dv = 0.5*(a(x)+a(x+dx))*dt
    return dx1,dx2,dv1,dv2 # 변화량 출력
    
winsize =500 #창 사이즈 지정
Ratom =0.5 # 원자 반지름
Matom = 1. # 원자 질량
epsilon = 10.2
sigma = 1.
l = 1.9
r_c = 2.5*sigma
Time = 0 # 초기 시간
dt = 0.001                                   # 시간 간격
x_v = np.arange(100.,500.,10.)
y_p = []
for y in x_v:
    atom = [] # 원자
    atom_pos = [] # 원자 위치
    atom_v = [] # 원자 속도
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
