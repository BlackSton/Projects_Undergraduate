import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg
from matplotlib import animation

def init():             # 그래프 초기화
    Wave.set_data([],[])
    Potential.set_data([],[])
    Time.set_text('')
    return(Wave,Potential,Time)
    
def animate(i):         # 그래프 그리기
    plotting = Result[:,i]
    Wave.set_data(x,plotting)
    Potential.set_data(x,np.abs(V))
    Time.set_text('time: '+ '%.2f'%(i*dt))
    return (Wave,Potential,Time)
    
def Calculation(J,N,psi0,V,dt,dx):  #수치 시뮬레이션 함수
    o = np.ones((J),complex)
    alp = (1j)*dt/(2*dx**2)*o
    xi = o + 1j*dt/2*(2/(dx**2)*o +V)
    gamma = o - 1j*dt/2*(2/(dx**2)*o +V)
    diags = np.array([-1,0,1])
    vecs1 = np.array([-alp,xi,-alp])
    vecs2 = np.array([alp,gamma,alp])
    U1 = spdiags(vecs1,diags,J,J)
    U2 = spdiags(vecs2,diags,J,J)
    U1 = U1.tocsc()
    U2 = U2.tocsc()
    PSI = np.zeros((J,N),complex)
    PSI[:,0] = psi0     #Put initial wave function
    LU = linalg.splu(U1)
    for n in range(0,N-1):
        b = U2.dot(PSI[:,n])
        PSI[:,n+1] = LU.solve(b)
    return PSI
# 기본 값 설정 구간
dx = 5E-2   # x 간격
dt = 1E-2   # 시간 간격
T  = 7      # 계산할 시간
x0 = -35    # 첫 x 좌표
x1 = 40     # 마지막 x 좌표
x = np.arange(x0,x1,dx) # 계산용 x값 생성
E = 5       #파동 에너지 값
k = np.sqrt(2*E)
V0 = 20
V = []
delta_x = 3.0
delta_E = k/(2*delta_x)

#함수 지정구간
psi0 = np.exp(-(x+10)**2/(2*delta_x**2))*np.exp(k*x*1j) #   초기 파동 함수
for i in x:                                             #   포텐셜함수 설정
    if i>5 and i<6:
        V.append(V0)
    else:
        V.append(0)

#시뮬레이션 구간
V = np.complex128(V) # 계산용 허수화
Result = Calculation(len(x),int(T/dt),psi0,V,dt,dx) #슈뢰딩거 방정식 풀기
Result = (V0/4)*np.abs(Result)**2 + V0/3

#그래프 그리는 구간
fig,ax = plt.subplots()          # 빈그래프 생성
ax.set_xlim((x0,x1))             # x범위 지정
ax.set_ylim((-V0*0.1,V0*1.1))    # y 범위 지정
Wave, = ax.plot([],[],lw =2,label="Probability")       #확률밀도함수 그래프
Potential, = ax.plot([],[],lw = 2,label = "Potential") #포텐셜함수 그래프
Time = ax.text(x1-15,V0, '')     #시간 출력용
ax.legend(loc=2)                 #범레 보여주기
anim = animation.FuncAnimation(fig,animate,init_func=init,
                               interval = 20, blit = True) #애니메이션화
plt.show()


