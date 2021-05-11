import numpy as np 
import matplotlib.pyplot as plt 
#
#
#
# Eigenvalues in x space 
data = np.loadtxt('eva_x.dat')
idxs = data[:,0]
#n_s = len(idxs)
n_s = 10
eva_x = data[:,1]

# Eigenfunctions in x_space 
data = np.loadtxt('eve_x.dat')
x = data[:,1]
V_x = data[:,2]
waves = data[:,3:]
fig, ax = plt.subplots()
#
_ = ax.set_title('Direct space eigenfunctions')
_ = ax.set_xlabel(r'$x$')
_ = ax.set_ylabel(r'$E$')
y_inf = np.amin(V_x) - 0.02 ; y_sup = eva_x[n_s -1] + 5.0 
_ = ax.set_ylim(y_inf, y_sup)
_ = ax.plot(x, V_x)
ff = ( eva_x[1] - eva_x[0] ) * 0.3
for i in range(n_s):
    _ = ax.axhline(y = eva_x[i], lw=0.5, color='cadetblue')
    _ = ax.plot(x, eva_x[i] + ff*waves[:,i] , color= 'cadetblue')
fig.savefig('waves_x.png')
#
#
#
# Eigenvalues comparison 
plt.close('all')
fig, ax = plt.subplots(1,2, figsize=(9,6))
plt.subplots_adjust(wspace=0.3)
data = np.loadtxt('eva_k.dat')
eva_k = data[:,1]
_ = ax[0].set_title(' Eigenvalues confrontation')
_ = ax[0].set_xlabel(' State index')
_ = ax[0].set_ylabel(r'$E/E_0$')
_ = ax[0].plot(idxs[:n_s], eva_x[:n_s], marker='.', color='cadetblue', label='Direct Eigenvalues', alpha=0.5)
_ = ax[0].plot(idxs[:n_s], eva_k[:n_s], marker='.', color='firebrick', label='Reciprocal eigenvalues', alpha=0.5)
_ = ax[0].legend()
_ = ax[1].set_title(' Absolute difference')
_ = ax[1].set_ylabel(r'$E_k - E_x$')
_ = ax[1].set_xlabel('State index')
_ = ax[1].plot(idxs[:n_s], (eva_k[:n_s] - eva_x[:n_s] ) , marker='.', color='coral', label='Difference')
fig.savefig('w_eva_confrontation.png')
#
#
#
# Fourier transform comparison
plt.close('all')
fig, ax = plt.subplots(3,2, figsize=(9,12))
plt.subplots_adjust(wspace=0.3)
data = np.loadtxt('four.dat')
x = data[:,0]
re_f = data[:,1]
im_f = data[:,2]
k = data[:,3]
re_F = data[:,4]
im_F = data[:,5]
re_at_f = data[:,6]
im_at_f = data[:,7]
#
_ = ax[0,0].set_ylabel(r'Re$(f)$')
_ = ax[0,0].set_xlabel(r'$x$')
_ = ax[0,1].set_ylabel(r'Im$(f)$')
_ = ax[0,1].set_xlabel(r'x')
_ = ax[0,0].plot(x, re_f, marker='.', color='cadetblue')
_ = ax[0,1].plot(x, im_f, marker='.', color='cadetblue')
#
_ = ax[1,0].set_ylabel(r'Re$(F)$')
_ = ax[1,0].set_xlabel(r'$k$')
_ = ax[1,1].set_ylabel(r'Im$(F)$')
_ = ax[1,1].set_xlabel(r'k')
_ = ax[1,0].plot(k[:n_s], re_F[:n_s], marker='.', color='firebrick')
_ = ax[1,1].plot(k[:n_s], im_F[:n_s], marker='.', color='firebrick')
#
_ = ax[2,0].set_ylabel(r'Re$(F^{-1})$')
_ = ax[2,0].set_xlabel(r'$x$')
_ = ax[2,1].set_ylabel(r'Im$(F^{-1})$')
_ = ax[2,1].set_xlabel(r'x')
_ = ax[2,0].plot(x, re_at_f, marker='.', color='cadetblue')
_ = ax[2,1].plot(x, im_at_f, marker='.', color='cadetblue')
fig.savefig('w_four.png')
#
#
#

plt.close('all')
data = np.loadtxt('eve_test.dat')
fig, ax = plt.subplots()
rebuild_f = data[:,0]
x_f = data[:,1]
_ = ax.set_xlabel(r'$x$')
_ = ax.set_ylabel(r'$|\psi_0(x)|^2$')
_ = ax.plot(x, rebuild_f, color='firebrick', label='Rebuilt eigenfunction')
_ = ax.plot(x, x_f, color='cadetblue', label='x space eigenfunction', alpha=0.75)
_ = ax.legend()
fig.savefig('w_test.png')






