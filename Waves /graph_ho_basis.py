import numpy as np 
import matplotlib.pyplot as plt 

#importing data from potentials.dat
data = np.loadtxt('potentials.dat')
x = data[:,0]
V_morse = data[:,1]
V_harmonic = data[:,2]

data = np.loadtxt('eva_mh.dat')
eva_m = data[:,0]
eva_h = data[:,1]
M = len(eva_m)

eve_morse = np.loadtxt('eve_morse.dat')
eve_harmonic = np.loadtxt('eve_harmonic.dat')

data = np.loadtxt('eva_eve_ho.dat')
eva_ho = data[:,0]
eve_ho = data[:,1:]

#potentials graph
fig, ax = plt.subplots(figsize=(9,6))
plt.subplots_adjust(wspace=0.4)
y_min = np.amin(V_morse); y_max = 0.5
_ = ax.set_ylim(y_min, y_max)
_ = ax.set_title('Potentials comparison')
_ = ax.set_xlabel(r'$x$')
_ = ax.set_ylabel(r'$V(x)$')
_ = ax.plot(x, V_morse, color='steelblue', label = 'Morse Eigenstates')
_ = ax.plot(x, V_harmonic, color = 'green', label= 'Harmonic Eigenstates')
for i in range(9) :
    _ = ax.axhline(y = eva_m[i], color='steelblue', alpha=0.5)
    _ = ax.plot(x, eva_m[i] + eve_morse[:,i], color='cadetblue', alpha=0.5)
    _ = ax.axhline(y = eva_h[i], color='green', alpha=0.5)
    _ = ax.plot(x, eva_h[i] + eve_harmonic[:,i], color='green', alpha=0.5)
_ = ax.legend()
fig.savefig('potentials.png')
#
#

fig, ax = plt.subplots(figsize=(9,6))
plt.subplots_adjust(wspace=0.4)
_ = ax.set_ylim(y_min, y_max)
_ = ax.set_xlim(-0.5, 10.5)
_ = ax.set_title('Eigenvalues conmparison')
_ = ax.set_xlabel('State index')
_ = ax.set_ylabel(r'$E (a. u.)$')
_ = ax.plot(range(M), eva_m, marker='.', color='steelblue', label = 'x space', alpha=0.5)
_ = ax.plot(range(M), eva_h, marker='.', color='green', label = 'x space but harmonic ones', alpha=0.5)
_ = ax.plot(range(M), eva_ho, marker='.', color='coral', label= 'Harmonic Oscillator space', alpha=0.5)
_ = ax.legend()
fig.savefig('eigenvalue comparison.png')
