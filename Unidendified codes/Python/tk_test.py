import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg
from matplotlib import animation

def init():             # reset graph
    line.set_data([],[])
    line1.set_data([],[])
    return(line,line1)
    
def animate(i):
    plotting = Result[:,i]
    line.set_data(x,plotting)
    line1.set_data(x,V)
    return (line,line1)
def make_zeta(V):
    m = 1
    q1 = 2*k
    q2 = 2*k
    q3 = 2*k
    h1 = m*(q1+q2+q3)
    h2 = (m**2)*(q1*q2*q3)*((1/q1)+(1/q2)+(1/q3))
    h3 = (m**3)*(q1*q2*q3)
    zeta1 = (-1j*(h2-V)/(2*dx)+1/(dt*dx)-1j*h1/(2*dt)-(h3-h2*V)/4)
    zeta2 = (1j*(h2-V)/(2*dx)-1/(dt*dx)-1j*h1/(2*dt)+(h3-h2*V)/4)
    zeta3 = (1j*(h2-V)/(2*dx)+1/(dt*dx)-1j*h1/(2*dt)+(h3-h2*V)/4)
    zeta4 = (-1j*(h2-V)/(2*dx)-1/(dt*dx)-1j*h1/(2*dt)+(h3-h2*V)/4)
    return np.array([0,zeta1,zeta2,zeta3,zeta4])
    
def Calculation(J,N,psi0,V):  #calculate schrodinger equation
    zeta = make_zeta(V[0])
    o = np.ones((J),complex)
    alp = (1j)*dt/(2*dx**2)*o
    xi = o + 1j*dt/2*(2/(dx**2)*o +V)
    gamma = o - 1j*dt/2*(2/(dx**2)*o +V)
    diags = np.array([-1,0,1])
    xi[0],xi[J-1] = zeta[1],zeta[1]
    gamma[0],gamma[J-1] = zeta[3],zeta[3]
    up1,up2,dn1,dn2 = -alp,alp,-alp,alp
    up1[1],dn1[J-2]=zeta[2],zeta[2]
    up2[1],dn2[J-2]=zeta[4],zeta[4]
    vecs1 = np.array([dn1,xi,up1])
    vecs2 = np.array([dn2,gamma,up2])
    U1 = spdiags(vecs1,diags,J,J)
    U2 = spdiags(vecs2,diags,J,J)
    U1 = U1.tocsc()
    U2 = U2.tocsc()
    PSI = np.zeros((J,N),complex)
    PSI[:,0] = psi0
    LU = linalg.splu(U1)
    for n in range(0,N-1):
        b = U2.dot(PSI[:,n])
        PSI[:,n+1] = LU.solve(b)
    return PSI
n= 2000
N = 400
dt = (10.-0.)/N
x = np.linspace(0,40,n)
dx = x[1]-x[0]
E = 5
k = np.sqrt(2*E)
V0 = -20
V = []
delta_x = 2.0
delta_E = k/(2*delta_x)
psi0 = np.exp(-(x-10)**2/(2*delta_x**2))*np.exp(k*x*1j) #   initial wave function
for i in x:                                             #   definition of Potential
    if i>15 and i<18:
        V.append(V0)
    else:
        V.append(0)
V = np.complex128(V) # to compucate wave function
fig,ax = plt.subplots() # create black plot
ax.set_xlim((0,40))   # set x limitation
ax.set_ylim((-1,10))    # set y limitaiton
line, = ax.plot([],[],lw = 2)
line1, = ax.plot([],[],lw = 2)
Result = Calculation(n,N,psi0,V)
Result = 2*np.abs(Result)**2 + 2
anim = animation.FuncAnimation(fig,animate,init_func=init,interval = 20, blit = True)
plt.show()
















