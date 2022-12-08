from visual import *
from visual.graph import *
from random import random

winsize =500
# 물리 상수 지정
L = 10
g = 0
Nball = 100
Ratom = 0.3
Matom = 1.
T = 300.
k = 1.4E-23
dt = 0.01

# 그래프 계산용 변수 지정
v=0
theory_sum = 0
count = 0

scene = display(title="PV", width=winsize, height=winsize, x=0, y=0,
                center=(L/2.,L/2.,L/2.))
vdist = gdisplay(x=0, y=winsize, ymax = 20,
             width=winsize, height=0.6*winsize, xtitle='v', ytitle='dN')
theory = gcurve(color=color.red)
theory_avg = gcurve(color = color.cyan)
graph = 0

floor_Top = box(pos=(L/2,L,L/2),size = (L,0.2,L),color = color.blue)
xaxis = curve(pos=[(0,0,0), (L,0,0)], color=color.blue)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=color.blue)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=color.blue)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=color.blue)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=color.blue)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=color.blue)

ball = []
ball_Pos = []
ball_Velocity = []
ball_p = []
for i in range(Nball):
    x = random()*L
    y = random()*L
    z = random()*L
    phi = random()*2*pi
    psi = random()*2*pi
    ball.append(sphere(pos=(x,y,z), radius=Ratom, color=color.red))
    pavg =sqrt(2.*Matom*1.5*k*T)
    ball_Pos.append((x,y,z))
    ball_p.append((pavg*sin(psi)*cos(phi),pavg*sin(psi)*cos(phi),pavg*cos(psi)))
Pos = array(ball_Pos)
p = array(ball_p)

floor_Top.velocity = vector(0,0,0)
while 1:
    v+= dt
    rate (100)
    Pos = Pos + (p/Matom)*dt
    floor_Top.velocity = floor_Top.velocity + vector(0,g*dt,0)
    floor_Top.pos = floor_Top.pos + floor_Top.velocity*dt
    outside = less_equal(Pos,Ratom) # walls closest to origin
    p1 = p*outside
    p = p-p1+abs(p1) # force p component inward
    graph = (sum(abs(2*p*outside)))
    if v < 5 and 0.5 < v:
        theory_sum += graph
        count += 1
        theory.plot(pos=(v-0.5,graph))
        theory_avg.plot(pos=(v-0.5,theory_sum/count))
        print(theory_sum/(count*3*L**2))
    outside = greater_equal(Pos,L-Ratom) # walls farther from origin
    p1 = p*outside
    p = p-p1-abs(p1) # force p component inward
    for i in range(Nball):
        ball[i].pos = Pos[i]

