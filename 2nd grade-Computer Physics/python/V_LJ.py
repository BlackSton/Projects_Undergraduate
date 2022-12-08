from __future__ import print_function,division
from visual import *
from visual.graph import *
winsize =500 #창 사이즈 지정
L = 10
Ratom =0.5
Matom = 1
epsilon = 0.01
sigma = 1.5
r_c = 2.5*sigma
Time = 0
dt = 0.0003
scene = display(title="12-6 potential", width=winsize, height=winsize,\
                x=0, y=0,center=(0.,0.,L/2.))
atom = [] # 원자
atom_pos = [] # 원자 위치
atom_v = [] # 원자 속도
atom_pos.append((10,0,0))
atom_pos.append((-10,0,0))
atom.append(sphere(pos=atom_pos[0],radius=Ratom,color = color.red))
atom.append(sphere(pos=atom_pos[1],radius=Ratom,color = color.red))
atom_v.append((-80,0,0))
atom_v.append((80,0,0))
pos = array(atom_pos)
v = array(atom_v)
a = 0
while True:
    Time += dt
    rate(1000)
    v[0][0] = v[0][0] + a*dt
    v[1][0] = v[1][0] - a*dt
    pos = pos + v*dt
    r = sqrt((pos[0][0]-pos[1][0])**2+(pos[0][1]-pos[1][1])**2) # 원자간의 거리
    if r < r_c:
        a = ((24*epsilon/sigma)*(2*(sigma/r)**13-(sigma/r)**7)-(24*epsilon/sigma)*(2*(sigma/r_c)**13-(sigma/r_c)**7))/Matom #
    else:
        a = 0
        if pos[0][0] > 10: 
            print(v[0][0],v[1][0])
            exit()
    atom[0].pos = pos[0]
    atom[1].pos = pos[1]
