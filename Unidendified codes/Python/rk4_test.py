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
    m = 0.5
    q1 = 2*k
    q2 = 2*k
    q3 = 2*k
    q4 = 2*k
    g1 = m*(q1+q2+q3+q4)
    g2 = (m**2)*((q1*q2)+(q1*q3)+(q1*q4)+(q2*q3)+(q2*q4)+(q3*q4))
    g3 = (m**3)*(q1*q2*q3*q4)*((1/q1)+(1/q2)+(1/q3)+(1/q4))
    g4 = (m**4)*q1*q2*q3*q4
    p1 = -1
    p2 = g1
    p3 = 1j*(g2-2*V)
    p4 = 1j*(g1*V-g3)
    p5 = V**2 - g2*V + g4
    zeta1 = p1/(2*dt**2)-p2/(dt*dx)+p3/(2*dt)-p4/(2*dx)+p5/4
    zeta2 = p1/(2*dt**2)+p2/(dt*dx)+p3/(2*dt)+p4/(2*dx)+p5/4
    zeta3 = p1/(2*dt**2)-p2/(dt*dx)+p3/(2*dt)+p4/(2*dx)-p5/4
    zeta4 = p1/(2*dt**2)+p2/(dt*dx)+p3/(2*dt)-p4/(2*dx)-p5/4
    zeta5 = -p1/(2*dt**2)
    zeta6 = -p2/(2*dt**2)
    return np.array([0,zeta1,zeta2,zeta3,zeta4,zeta5,zeta6])
    
def Calculation(J,N,psi0,V):  #calculate schrodinger equation
    zeta = make_zeta(V)
    o = np.ones((J),complex)
    alp = (1j)*dt/(2*dx**2)*o
    xi = o + 1j*dt/2*(2/(dx**2)*o +V)
    gamma = o - 1j*dt/2*(2/(dx**2)*o +V)
    diags = np.array([-1,0,1])
    vecs1_1 = np.array([-alp,xi,-alp])
    vecs2_1 = np.array([alp,gamma,alp])
    W1 = spdiags(vecs1_1,diags,J,J)
    W2 = spdiags(vecs2_1,diags,J,J)
    W1.tocsc()
    W2.tocsc()
    xi[0],xi[J-1] = zeta[1],zeta[1]
    gamma[0],gamma[J-1] = zeta[3],zeta[3]
    up1,up2,dn1,dn2 = -alp,alp,-alp,alp
    up1[1],dn1[J-2]=zeta[2],zeta[2]
    up2[1],dn2[J-2]=zeta[4],zeta[4]
    vecs1 = np.array([dn1,xi,up1])
    vecs2 = np.array([dn2,gamma,up2])
    U1 = spdiags(vecs1_1,diags,J,J)
    U2 = spdiags(vecs2_1,diags,J,J)
    U1 = U1.tocsc()
    U2 = U2.tocsc()
    Z  = np.zeros((J,J),complex)
    Z[0,0] = zeta[5]
    Z[0,1] = zeta[6]
    Z[-1,-1] = zeta[5]
    Z[-1,-2] = zeta[6]
    PSI = np.zeros((J,N),complex)
    PSI[:,0] = psi0
    LU = linalg.splu(U1)
    LW = linalg.splu(W1)
    for n in range(0,N-1):
        if n == 0:
            b = W2.dot(PSI[:,n])
            PSI[:,n+1] = LW.solve(b)
        else:
            b = U2.dot(PSI[:,n])# + Z.dot(PSI[:,n-1])
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
ax.set_xlim((10,30))   # set x limitation
ax.set_ylim((-1,6))    # set y limitaiton
line, = ax.plot([],[],lw = 2)
line1, = ax.plot([],[],lw = 2)
Result = Calculation(n,N,psi0,V)
Result = 2*np.abs(Result)**2 + 2
anim = animation.FuncAnimation(fig,animate,init_func=init,interval = 20, blit = True)
plt.show()
















