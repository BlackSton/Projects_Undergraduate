from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

winsize =500 #창 사이즈 지정

# 일반 상수 지정
L = 10 # 용기 한변 길이(단위:m)
B_L = 0.5 # 윗벽 두깨(단위:m)
Natom = 30 #원자 갯수
Ratom = 0.1  #원자 반지름(단위:m)
Matom = 4E-3/6E23 # 원자 질량 4g/아보가드로수 = 0.004/6.02*10^23
g = 9.8 # 중력상수
k = 1.4E-23 # 볼츠만 상수
dt = 1.E-4 #시간 간격
beta =1. # 마찰계수
T =300 #온도(300K)

# 그래프 계산용 변수 지정
time=0 # 초기 시간
theory_sum = 0 #평균 압력측정용
theory_avg = [] #평균 압력-그래프-
M_B_Min = 1E-20 #최소 무게
M_B_Max = 5E-20 #최고 무게
dM = 5E-21 #무게 간격
graph_color =[] #그래프 색깔

## 그래프 설정 란 ##
Time_to_P = gdisplay(title="시간에 따른 압력",x = 0,y = winsize,ymax = M_B_Max/500,\
                  width = winsize, height = 0.6*winsize, xtitle ="Time",ytitle="P") #시간에 따른 압력 그래프 화면 생성
Time_P = gdots(color=color.red)#시간에 따른 압력 그래프 점 찍는거
T_to_P = gdisplay(title="압력과 부피간의 관계",x = winsize,y = 0,ymax = L**3,\
                  width = winsize, height = winsize, xtitle ="P",ytitle="V") #온도와 압력관게  그래프화면 생성
P_T_Theory = gcurve(color = color.cyan) #이론 그래프 색깔 지정
P_T_Experimental = gdots(color=color.red) #측정 그래프 색깔 지정
for graph_M in arange(dM,M_B_Max+1E-22+dM,dM): # 이론 예측 그래프 생성
    P_T_Theory.plot(pos=(graph_M,Natom*k*L**3/(graph_M*9.8)))

time_to_P = gdisplay(title="시간에 따른 부피",x=winsize, y=winsize,\
                     ymax = L**3,width=winsize, height=0.6*winsize, xtitle='Time', ytitle='V') #시간에 따른 부피 그래프 화면 생성
theory = gdots(color=color.red)
for i in range(int((M_B_Max-M_B_Min+dM)/dM)):
    graph_color.append((0.5,i/8.,i/8.+0.3))
    theory_avg.append(gcurve(color =graph_color[i]))

for M in arange(M_B_Min,M_B_Max+1E-22,dM):  #M_B_Min 에서 B_M_Max 까지 dM 간격으로 시뮬레이션..
    ## 시뮬레이션 창생성 ##
    scene = display(title="압력에 따른 부피", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

    ## 외벽 선 생성 ##
    xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.blue)
    yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.blue)
    zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.blue)
    xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.blue)
    yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.blue)
    zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.blue)

    ## 천장 생성 ##
    T_B = box(pos=vector(L/2,L,L/2),size=vector(L,B_L,L))
    B_Pos = vector(L/2,L,L/2)
    B_v = vector(0,0,0)
    
    atom = [] #원자 
    atom_Pos = [] #원자 위치
    atom_p = [] # 원자 운동량

    ## 원자의 위치와 속도 지정##
    for i in range(Natom):
        ## 원자의 위치 랜덤화 ##
        x = Ratom+random()*(L-2*Ratom)
        y = Ratom+random()*(L-2*Ratom)
        z = Ratom+random()*(L-2*Ratom)
        atom_Pos.append((x,y,z))
        atom.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))

        ## 속도값 지정과 방향 랜덤화 ##
        phi = random()*2*pi
        psi = random()*2*pi
        pavg =sqrt(3.*Matom*k*T) #원자의 평균 속도 지정
        atom_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*sin(phi),pavg*cos(psi)))

    Pos = array(atom_Pos)
    p = array(atom_p)
    theory_length = 500 # 평균할 범위 지정
    theory_P_avg = empty([theory_length],float)# 평균을위한 값 생성
    theory_P_avg[:]=1E-23# 평균 초기값(온도에따라 설정해줘야 함.<-첫 평균압력값임 그래프 보았을때 이상하면 이값 수정하면 됨.)
    while time < 3.0:
        time+= dt
        rate (1000)
        ## dt의 시간이 지난후 위치 계산 ##
        Pos = Pos + (p/Matom)*dt
        B_v = B_v - vector(0,(9.8)*dt,0)
        ## 윗 면들 충돌 계산 ##
        outside = greater_equal(Pos,[L-Ratom,B_Pos[1]-Ratom-B_L,L-Ratom])
        p1 = p*outside
        p = p-p1-abs(p1)
        B_V = B_v[1] #충돌전 벽의 속도
        theory_P = 0
        for c in range(Natom):
            if outside[c][1]:
                b_v = -p[c][1]/Matom # 충돌전 공의 속도
                theory_P += b_v*Matom
                b_v0= (b_v*(Matom-M)+2*M*B_V)/(Matom+M)#충돌후 공의 속도
                B_V0= (B_V*(M-Matom)+2*Matom*b_v)/(Matom+M)#충돌후 벽의 속도
                B_V = B_V0
                theory_P += b_v *Matom
                p[c][1] = b_v0*Matom
        B_v = vector(0,B_V,0)
        B_Pos = B_Pos + B_v*dt
        ## 시간에 따른 압력 값 출력(윗면의 압력만 출력하는것임)##
        theory_P_avg=append(theory_P_avg,theory_P)
        theory_P_avg=delete(theory_P_avg,0)
        Time_P.plot(pos=(time,sum(theory_P_avg)/theory_length),color=graph_color[int((M-M_B_Min)/dM)])
        ## 밑바닥 면들 충돌 계산 ##
        outside = less_equal(Pos,Ratom)
        p1 = p*outside
        p = p-p1+abs(p1)
        ## 시간에 따른 부피 값 출력 ##
        theory_volume= L*L*(B_Pos[1]-B_L)
        theory_avg[int((M-M_B_Min)/dM)].plot(pos=(time,theory_volume))
        ## 원자 위치 재설정 ##   
        for i in range(Natom):
            atom[i].pos = Pos[i]
        T_B.pos = B_Pos
    P_T_Experimental.plot(pos=(M,theory_volume),color=graph_color[int((M-M_B_Min)/dM)]) #측정값 그래프화

    ## 값 초기화 ##
    time = 0
    theory_sum = 0
    scene.delete()

while 1:
    rate(1)
