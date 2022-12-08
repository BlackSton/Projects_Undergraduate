import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags,linalg
from matplotlib import animation


def init():             # reset graph
    line.set_data([],[])
    return(line,)
    
def animate(i):
    plotting = Result[:,i]
    line.set_data(x,plotting)
    return (line,)
    
def Calculation(J,N,psi0,V,dt,dx):  #calculate schrodinger equation
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
    PSI[:,0] = psi0
    LU = linalg.splu(U1)
    for n in range(0,N-1):
        b = U2.dot(PSI[:,n])
        PSI[:,n+1] = LU.solve(b)
    return PSI
n = 4000
N = 500
dt = (7.-0.)/N
x = np.linspace(-35,40,n)
dx = x[1]-x[0]
E = 5
k = np.sqrt(2*E)
V0 = 10
V = []
delta_x = 2.0
delta_E = k/(2*delta_x)
psi0 = np.exp(-(x+10)**2/(2*delta_x**2))*np.exp(k*x*1j) #   initial wave function
for i in x:                                             #   definition of Potential
    if i>5 and i<8:
        V.append(V0)
    else:
        V.append(0)
V = np.complex128(V) # to compucate wave function
fig,ax = plt.subplots() # create black plot
ax.set_xlim((-35,40))   # set x limitation
ax.set_ylim((-5,5))    # set y limitaiton
line, = ax.plot([],[],lw =2)

Result = Calculation(n,N,psi0,V,dt,dx)
anim = animation.FuncAnimation(fig,animate,init_func=init,frames = np.arange(0,N),interval = 100, blit = True)
anim.save('Wave_.gif', writer='imagemagick', fps=30)







