import numpy as np

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

with open('../Input/data_sheet.dat', 'r') as f:
    for i in range(2):
        f.readline()
    L = f.readline().split()
    L = float(L[1])
    for i in range(16):
        f.readline()
    save_wave =f.readline().split()
    save_wave = int(save_wave[1]) + 1


wave_func = np.loadtxt('../Output/moving.txt')
V_f = np.loadtxt('../Output/pot.txt')
V = V_f[:,1]

x = wave_func[:,0]

fig = plt.figure()
ax = plt.axes(xlim=(0,L), ylim=(-1,1))
ax.plot(x, V)
ln, = ax.plot([], [])

def ani(t):
    ln.set_data(x, wave_func[:,2*t+1])
    return ln,

animat = FuncAnimation(fig, ani, frames=save_wave, interval=20, blit=True)
plt.show()



