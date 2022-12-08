import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()
ax.set_xlim((0,2))
ax.set_ylim((-2,2))
ax.grid(True)

line, = ax.plot([],[],lw=2)
def init():
    line.set_data(([],[]))
    return (line, )

def animate(t):
    x = np.linspace(0,2,1000)
    y = np.sin(2*np.pi*(x-0.01*t))
    line.set_data(x,y)
    print(t)
    return (line, )

ani = animation.FuncAnimation(fig=fig, func=animate, init_func = init, interval = 1000,blit = True)
plt.show()