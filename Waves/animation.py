import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

with open('Input/data_sheet.dat', 'r') as f:
    for i in range(2):
        f.readline()
    L = f.readline().split()
    L = float(L[1])
    for i in range(23):
        f.readline()
    save_wave = f.readline().split()
    save_wave = int(save_wave[1]) + 1


wave_func = np.loadtxt('Output/moving.txt')
V_f = np.loadtxt('Output/pot.txt')
x = V_f[:,0]
V = V_f[:,1]

fig = plt.figure()
ax = plt.axes(xlim=(0,L), ylim=(0.,np.max(wave_func[:, 1:])))
ax.plot(x, V)
ln, = ax.plot([], [])

def ani(t):
    ln.set_data(x, wave_func[:,t])
    return ln,

animat = FuncAnimation(fig, ani, frames=save_wave, interval=20, blit=True)
plt.show()



