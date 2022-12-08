from __future__ import print_function,division
from visual import *
from visual.graph import *
from math import sqrt
import numpy as np

def distance(pos1,pos2):
    return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)

def V_LJ_a(pos1,pos2):
    global epsilon,r_c,sigma,Matom
    a = np.array([0.,0.,0.])
    r_d = distance(pos1,pos2)
    if r_d < r_c:
        r = pos1 - pos2
        a = ((24*epsilon*r/sigma**2)*(2*(sigma/r_d)**14-(sigma/r_d)**8)\
             -(24*epsilon*r/sigma**2)*(2*(sigma/r_c)**14-(sigma/r_c)**8))/Matom
    else:
        a = np.array([0.,0.,0.])
    return a
def Verlet(pos1,pos2,v1,v2,a):
    global dt
    dx1 = (v1 * dt) + (0.5 * a(pos1,pos2) * dt**2)
    dx2 = (v2 * dt) + (0.5 * a(pos2,pos1) * dt**2)
    V_a_dx1 = a(pos1+dx1,pos2+dx2)
    V_a_dx2 = a(pos2+dx2,pos1+dx1)
    dv1 = 0.5*(a(pos1,pos2) + V_a_dx1)*dt
    dv2 = 0.5*(a(pos2,pos1) + V_a_dx2)*dt
    return dx1,dx2,dv1,dv2

winsize =500 #창 사이즈 지정

Natoms = 100  # change this to have more or fewer atoms

# Typical values
L = 1. # container is a cube L on a side
gray = (0.7,0.7,0.7) # color of edges of container
Matom = 4E-3/6E23 # helium mass
Ratom = 0.03 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant
T = 300. # around room temperature
dt = 1E-5

epsilon = 10.2
sigma = 1.
r_c = 2.5*sigma

scene = display(title="12-6 potential", width=winsize, height=winsize,\
                x=0, y=0,center=(L/2.,L/2.,L/2.))
deltav = 100. # binning for v histogram
vdist = gdisplay(x=0, y=winsize, ymax = Natoms*deltav/1000.,
             width=winsize, height=0.6*winsize, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan)

dv = 10.
for v in arange(0.,3001.+dv,dv): # theoretical prediction
    theory.plot(pos=(v,
        (deltav/dv)*Natoms*4.*pi*((Matom/(2.*pi*k*T))**1.5)
                     *exp((-0.5*Matom*v**2)/(k*T))*v**2*dv))


xaxis = curve(pos=[(0,0,0), (L,0,0)], color=gray)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=gray)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=gray)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=gray)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=gray)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=gray)

Atoms = []
##colors = [color.red, color.green, color.blue,
##          color.yellow, color.cyan, color.magenta]
poslist = []
plist = []
mlist = []
rlist = []

for i in range(Natoms):
    Lmin = 1.1*Ratom
    Lmax = L-Lmin
    x = Lmin+(Lmax-Lmin)*random()
    y = Lmin+(Lmax-Lmin)*random()
    z = Lmin+(Lmax-Lmin)*random()
    r = Ratom
    if i == 0:
        Atoms.append(sphere(pos=(x,y,z), radius=r, color=color.yellow,
                    make_trail=True, retain=100))
    else:
        Atoms.append(sphere(pos=(x,y,z), radius=r, color=color.cyan))
    mass = Matom*r**3/Ratom**3
    pavg = sqrt(2.*mass*1.5*k*T) # average kinetic energy p**2/(2mass) = (3/2)kT
    theta = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    poslist.append((x,y,z))
    plist.append((px,py,pz))
    mlist.append(mass)
    rlist.append(r)

pos = array(poslist)
p = array(plist)
m = array(mlist)
m.shape = (Natoms,1) # Numeric Python: (1 by Natoms) vs. (Natoms by 1)
radius = array(rlist)

while True:
    rate(50)
    observation.plot(data=mag(p/m))
    # Update all positions
    r = pos-pos[:,newaxis] # all pairs of atom-to-atom vectors
    rmag = sqrt(sum(square(r),-1)) # atom-to-atom scalar distances
    hit = less_equal(rmag,radius+radius[:,None])-identity(Natoms)
    hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*Natoms+j

    # If any collisions took place:
    for ij in hitlist:
        i, j = divmod(ij,Natoms) # decode atom pair
        hitlist.remove(j*Natoms+i) # remove symmetric j,i pair from list
        dx0,dx1,dv0,dv1 = Verlet(pos[i],pos[j],v[i],v[j],V_LJ_a)
        v[i] = v[i] + dv0
        v[j] = v[j] + dv1
        pos[i] = pos[i] + dx0
        pos[j] = pos[j] + dx1        
 
    # Bounce off walls
    outside = less_equal(pos,Ratom) # walls closest to origin
    p1 = p*outside
    p = p-p1+abs(p1) # force p component inward
    outside = greater_equal(pos,L-Ratom) # walls farther from origin
    p1 = p*outside
    p = p-p1-abs(p1) # force p component inward

    # Update positions of display objects
    for i in range(Natoms):
        Atoms[i].pos = pos[i]

