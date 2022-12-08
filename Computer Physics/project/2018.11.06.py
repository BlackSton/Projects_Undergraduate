from __future__ import print_function,division
from visual import *
from visual.graph import *
from random import random

winsize =500 #창 사이즈 지정

# 일반 상수 지정
L = 10 # 용기 한변 길이
Natom = 30 #원자 갯수
Ratom = 0.4  #원자 반지
Matom = 4E-3/6E23 # 원자 질량 4g/아보가드로수 = 0.004/6.02*10^23
k = 1.4E-23 # 볼츠만 상수
dt = 1.E-4 #시간 간격

# 그래프 계산용 변수 지정
time=0 # 초기 시간
theory_sum = 0 #평균 압력측정용
theory_avg = [] #평균 압력-그래프-
T_Min = 50 #최소 온도
T_Max = 500 #최고 온도
dT = 50. #온도 간격
graph_color =[] #그래프 색깔

## 그래프 설정 란 ##
T_to_P = gdisplay(title="온도와 압력간의 관계",x = winsize,y = 0,ymax = T_Max*Natom*k/L**3,\
                  width = winsize, height = winsize, xtitle ="T",ytitle="P") #온도와 압력관게  그래프화면 생성
P_T_Theory = gcurve(color = color.cyan) #이론 그래프 색깔 지정
P_T_Experimental = gdots(color=color.red) #측정 그래프 색깔 지정
for graph_T in arange(0.,501.+dT,dT): # 이론 예측 그래프 생성
    P_T_Theory.plot(pos=(graph_T,graph_T*Natom*k/L**3))

time_to_P = gdisplay(title="시간에 따른 평균 압력",x=winsize, y=winsize,\
                     ymax = T_Max*Natom*k/L**3,width=winsize, height=0.6*winsize, xtitle='Time', ytitle='p') #시간에 따른 평균 압력 그래프화면 생성
theory = gdots(color=color.red)
for i in range(int((T_Max-T_Min+dT)/dT)):
    graph_color.append((0.5,i/8.,i/8.+0.3))
    theory_avg.append(gcurve(color =graph_color[i]))

for T in arange(T_Min,T_Max+1,dT):  #T_Min 에서 T_Max 까지 dT 간격으로 시뮬레이션..

    ## 시뮬레이션 창생성 ##
    scene = display(title="온도에 따른 압력", width=winsize, height=winsize, x=0, y=0,center=(L/2.,L/2.,L/2.))

    ## 외벽 선 생성 ##
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

    while time < 0.04:
        time+= dt
        rate (100)
        ## dt의 시간이 지난후 위치 계산 ##
        Pos = Pos + (p/Matom)*dt

        ## 밑바닥 면들 충돌 계산 ##
        outside = less_equal(Pos,Ratom)
        p1 = p*outside
        p = p-p1+abs(p1)

        ## 밑바닥 충돌시 압력 계산과 그래프 출력 ##
        theory_sum += (sum(abs(2*p*outside)))/dt
        theory_sum_avg = theory_sum/((time/dt)*3*L**2)
        #theory.plot(pos=(time,graph/(dt*3*L**2)))
        theory_avg[int((T-T_Min)/dT)].plot(pos=(time,theory_sum_avg))

        ## 윗 면들 충돌 계산 ##
        outside = greater_equal(Pos,L-Ratom)
        p1 = p*outside
        p = p-p1-abs(p1)
        
        ## 원자 위치 재설정 ##   
        for i in range(Natom):
            atom[i].pos = Pos[i]
   
    print("온도:",T,"측정 평균 압력:",theory_sum_avg,"이론 평균 압력:",T*Natom*k/L**3,"오차율",100*abs(theory_sum_avg-T*Natom*k/L**3)/(T*Natom*k/L**3),"%")

    P_T_Experimental.plot(pos=(T,theory_sum_avg),color=graph_color[int((T-T_Min)/dT)]) #측정값 그래프화

    ## 값 초기화 ##
    time = 0
    theory_sum = 0
    scene.delete()

while 1:
    rate(1)
