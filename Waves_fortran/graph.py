import numpy as np 
import matplotlib.pyplot as plt

data = np.loadtxt('Output/pot.txt')
x = data[:,0]
V_morse = data[:,1]

data = np.loadtxt('Output/pot_ho.txt')
V_ho = data[:,1] 

data = np.loadtxt('Output/eigenvalues_ho.txt')
indx = data[:,0]
eva_ho_space = data[:,1]

data = np.loadtxt('Output/eigenvalues.txt')
eva_x_space = data[:,1]


fig, ax = plt.subplots(1,2, figsize=(9,6))
plt.subplots_adjust(wspace=0.5)
ymin = np.amin(V_morse)-0.5; ymax=0.
_ = ax[0].set_ylim(ymin,ymax) 
_ = ax[0].plot(x, V_morse, label='Morse', color = 'cadetblue')
_ = ax[0].plot(x, V_ho, label='Harmonic approximation', color = 'coral')
_ = ax[0].set_xlabel(r'$r (a. u.)$')
_ = ax[0].set_ylabel(r'$V(r) (a. u.)$')
_ = ax[0].legend()

_ = ax[1].plot(indx, eva_x_space, label = 'R space eigenvalues', marker='.', color='cadetblue', alpha=0.5)
_ = ax[1].plot(indx, eva_ho_space, label= 'HO space eigenvalues', marker='.', color= 'coral', alpha=0.75)
_ = ax[1].set_xlabel('i (Eigenstate index)')
_ = ax[1].set_ylabel(r'$E_i (a.u)$')
_ = ax[1].legend()

fig.savefig('Images/test.png')