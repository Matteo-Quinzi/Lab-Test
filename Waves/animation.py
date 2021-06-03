import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

wave_func = np.loadtxt('Output/moving.txt')   # read the evolution file
num_frames = wave_func.shape[1]   # number of saved waves 

V_f = np.loadtxt('Output/pot.txt')   # read potential file
x = V_f[:,0]   # spatial interval
L = np.max(x)
V = V_f[:,1]   # potential
V0 = np.max(V)
V[:] = V[:] / V0 * 1.1 * np.max(wave_func[:,0])   # graphical purposes

fig = plt.figure()
ax = plt.axes(xlim=(0,L), ylim=(0.,np.max(V)*1.1))

ax.set_xlabel('x', fontsize='large')
ax.set_ylabel(r'$|\psi|^2$', fontsize='large')
secax = ax.secondary_yaxis('right')
secax.set_ylabel('Potential', fontsize='large')   # potential energy axis
secax.set_yticks(np.linspace(0., np.max(V), 4))
secax.set_yticklabels(np.around(np.linspace(0., V0, 4), decimals=2))

def ani(t):   # animation function
    ln.set_data(x, wave_func[:,t])
    return ln,

ax.plot(x, V)
ln, = ax.plot([], [])
animat = FuncAnimation(fig, ani, frames=num_frames, interval=20, blit=True)
plt.show()



