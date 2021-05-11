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
fig, ax = plt.subplots(1,2, figsize=(9,6))
plt.subplots_adjust(wspace=0.4)
y_min = np.amin(V_morse); y_max = 0.5
_ = ax[0].set_ylim(y_min, y_max)
_ = ax[0].set_title('Potentials comparison')
_ = ax[0].set_xlabel(r'$x$')
_ = ax[0].set_ylabel(r'$V(x)$')
_ = ax[0].plot(x, V_morse, color='steelblue', label = 'Morse Eigenstates')
_ = ax[0].plot(x, V_harmonic, color = 'green', label= 'Harmonic Eigenstates')
for i in range(6) :
    _ = ax[0].axhline(y = eva_m[i], color='steelblue', alpha=0.5)
    _ = ax[0].plot(x, eva_m[i] + eve_morse[:,i], color='cadetblue', alpha=0.5)
    _ = ax[0].axhline(y = eva_h[i], color='green', alpha=0.5)
    _ = ax[0].plot(x, eva_h[i] + eve_harmonic[:,i], color='green', alpha=0.5)
_ = ax[0].legend()
#
#
_ = ax[1].set_ylim(y_min, y_max)
_ = ax[1].set_xlim(-0.5, 10.5)
_ = ax[1].set_title('Eigenvalues conmparison')
_ = ax[1].set_xlabel('State index')
_ = ax[1].set_ylabel(r'$E (a. u.)$')
_ = ax[1].plot(range(M), eva_m, marker='.', color='steelblue', label = 'x space', alpha=0.5)
_ = ax[1].plot(range(M), eva_h, marker='.', color='green', label = 'x space but harmonic ones', alpha=0.5)
_ = ax[1].plot(range(M), eva_ho, marker='.', color='coral', label= 'Harmonic Oscillator space', alpha=0.5)
_ = ax[1].legend()
fig.savefig('potentials.png')
