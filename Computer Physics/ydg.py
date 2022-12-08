from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random
from math import acos,sqrt

def distance(pos1,pos2):
    return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)
def c1_select(c): # c1 값 반환 해주는 함수
    global Lc
    a_s = [] # c1값 배열
    for i in range(3):
        a_z = -Lc**2 + i*Lc**2 # Z축 +- 시키는것
        a = arange(a_z+c-Lc-1,a_z+c-Lc+2) # x축값 만듬
        a1 = append(a,a+Lc) # y축 +
        a_2 = append(a1,a+Lc*2) #y축 ++
        a_s = append(a_s,a_2) #z축에 추가
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

winsize =500 #창 사이즈 지정
epsilon = 10.2
sigma = 0.5
r_c = 2.5*sigma
# 일반 상수 지정
L = 10 # 용기 한변 길이(단위:m)
natom = 5 # 한 변당 원자갯수
Natom = natom**3 #총원자 갯수
Ratom = 0.5  #원자 반지름(단위:m)
Matom = 1. # 원자 질량 4g/아보가드로수 = 0.004/6.02*10^23
k = 1. # 볼츠만 상수
dt = 1.E-4 #시간 간격
T =100. #온도(300K)
## 시뮬레이션 창생성 ##
scene = display(title="압력에 따른 부피", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.blue)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.blue)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.blue)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.blue)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.blue)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.blue)

atom = [] #원자 
atom_Pos = [] #원자 위치
atom_p = [] # 원자 운동량

## 원자의 위치와 속도 지정##
## 원자의 위치 격자 구조 지정##
for i in range(natom):
    for j in range(natom):
        for K in range(natom):
            x = ((L-7*Ratom)/(natom-1))*i+3.5*Ratom # 벽과의 간격이 3.5*원자반지름으로 만든다.
            y = ((L-7*Ratom)/(natom-1))*j+3.5*Ratom # 벽과의 간격이 3.5*원자반지름으로 만든다.
            z = ((L-7*Ratom)/(natom-1))*K+3.5*Ratom # 벽과의 간격이 3.5*원자반지름으로 만든다.
            atom_Pos.append((x,y,z))
            atom.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))
            ## 속도값 지정과 방향 랜덤화 ##
            phi = random()*2*pi
            psi = acos(1-2*random())        # 보정값
            pavg = sqrt(2.*Matom*1.5*k*T) #원자의 평균 속도 지정
            atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

Pos = array(atom_Pos)       # 각 원자 배치 어레이
v = array(atom_p)/Matom     # 각 원자 속도 어레이


Lc=int(2+L//r_c)            # 한축에서의 칸의 갯수    +2 하는 이유는 벽 외벽라인에서도 계산해주기 위해서.
mc=empty(3,int)             # 한축에서의 칸 위치
lscl=zeros(Natom,int)           # lscl 칸만들기
head=empty(Lc**3,int) 

while True:
    rate(1000)
    #연결리스트 생성#
    for c in range(Lc**3):
        head[c]=-1
    for i in range(Natom):
        for a in range(3):
            mc[a]=floor(Pos[i,a]/(L/Lc))      # 모든 부분에서 설명 가능하게 하기위해바꿈
        c=mc[0]*Lc**2+mc[1]*Lc+mc[2]
        lscl[i]=head[c]                       # 뒤에 head 받아먹기
        head[c]=i
        
    # 상호 작용 하지 않는 원자들은 그냥 움직인다.
    for cc in range(Natom):
        if lscl[cc] == -1:
            Pos[cc] = Pos[cc] + v[cc]*dt
            
    #상호작용 계산 코드#
    for c in range(Lc**3):
        i = head[c]
        while i!= -1:
            for c1 in c1_select(c):
                j = head[int(c1)]
                while j != -1:
                    if i < j:
                        if distance(Pos[i],Pos[j])<r_c:
                            #상호작용
                            dx0,dx1,dv0,dv1 = Verlet(Pos[i],Pos[j],v[i],v[j],V_LJ_a)
                            v[i] = v[i] + dv0
                            v[j] = v[j] + dv1
                            Pos[i] = Pos[i] + dx0
                            Pos[j] = Pos[j] + dx1
                    j = lscl[j]
            i = lscl[i]
    
        ## 윗 면들 충돌 계산 ##
    outside = greater_equal(Pos,L-Ratom-L/Lc) # L/Lc는 한칸의 길이
    v1 = v*outside
    v = v-v1-abs(v1)
    ## 밑바닥 면들 충돌 계산 ##
    outside = less_equal(Pos,Ratom+L/Lc) # L/Lc는 한칸의 길이
    v1 = v*outside
    v = v-v1+abs(v1)
    for i in range(Natom):
        atom[i].pos = Pos[i]
