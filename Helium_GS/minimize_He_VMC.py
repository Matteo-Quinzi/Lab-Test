from HeVMC import *
import numpy as np 
import matplotlib.pyplot as plt
from time import process_time

#Reading data from 
input_file = 'He_VMC_input.txt'
fh = open(input_file,'r')
data = fh.readlines()

try :
    #Helium parameters
    interacting, jastrow, numeric = tuple(map(int,data[2].split()))
    try:
        interacting = bool(interacting)
        jastrow = bool(jastrow)
        numeric = bool(numeric)
    except:
        raise ValueError('Only 0 or 1 can be accepted as boolean parameters.')
    #metropolis parameters
    data_metropolis = data[6].split()
    steps = int(data_metropolis[0])
    d = float(data_metropolis[1])
    #statistical analysis parameters 
    data_statistics = data[10].split()
    method = int(data_statistics[0]) 
    if method == 0:
        method = 'split'
    elif method == 1: 
        method = 'bootstrap'
    else :
        raise ValueError("Invalid input for 'method' argument.")
    samples = int(data_statistics[1])
    #minimum implementation
    data_param = data[14].split()
    b_or_z = int(data_param[0]) 
    if b_or_z == 0 : 
        b_or_z = 'z'
    elif b_or_z == 1 :
        b_or_z = 'b'
    else:
        raise ValueError('Invalid input. Insert 0 or 1 into b_or_z argument.')
    param_min, param_max = tuple(map(float,data_param[1:3]))
    param_values = int(data_param[-1])
    #ds values 
    n_ds = int(data[18].split()[0])
    #gss values 
    gs_data = data[22].split()
    n_gs = int(gs_data[0])
    tol = float(gs_data[1])
    maxit = int(gs_data[2])

    #random number generator seed 
    seed = int(data[25].split()[0])
    
    #frequency to save steps into chain.txt
    freq = int(data[28].split()[0])

except: 
    raise ValueError('Something went wrong during input reading')

fh.close()

#Uncomment the next line to set the seed
np.random.seed(seed)


# A direct section search is implemented to 
# determine the region of the minimum 
#
if b_or_z == 'b': 
    b = np.linspace(param_min, param_max, param_values)
    Z = np.full(param_values,2.)
elif b_or_z == 'z':
    b = np.zeros(param_values)
    Z = np.linspace(param_min, param_max, param_values)

He_params = {
    'Z' : Z[0],
    'b' : b[0],
    'interacting' : interacting,
    'jastrow' : jastrow,
    'numeric' : numeric
}

energy_params = {
    'steps' : steps,
    'd' : d,
    'method' : method, 
    'samples' : samples
}

He = Helium(**He_params)
ds_energies = np.empty((n_ds,param_values)) 
ds_errors = np.empty((n_ds,param_values))
    
tic = process_time()
print('Performing Direct Search...\n')
for i in range(n_ds):
    _, _, ds_energies[i,:], ds_errors[i,:] = optimize_He_ds(Z, b, He, energy_args = energy_params)
    print('DS Run %d completed' % (i+1))
toc = process_time()

print('Direct search time : %7.5f s' % (toc-tic))

#Averaging quantities to reduce statistical noise 
avg_ds_energies = np.mean(ds_energies, axis=0)
avg_ds_errors = np.mean(ds_errors, axis=0)

#both informations on energy functional value and associated error are used 
#to shrink the interval before calling the golden section search method
#in order to ensure convergence in a reasunable interval
min_indx_en = np.argmin(avg_ds_energies)
min_indx_err = np.argmin(avg_ds_errors)
low_bound = min(min_indx_en, min_indx_err)
up_bound = max(min_indx_en,min_indx_err)
if low_bound != 0 : low_bound -= 1
if (up_bound  != len(b) - 1) :  up_bound += 1
    
if b_or_z == 'b' :
    param_bounds = (b[low_bound], b[up_bound])
elif b_or_z == 'z' :
    param_bounds = (Z[low_bound], Z[up_bound])

print('Starting Golden Section Search in (%5.4f,%5.4f)' % (param_bounds[0],param_bounds[1]) )
#Golden section search is called 
best_param = np.empty(n_gs)
print('Performing Golden Section Search ...')
tic = process_time()
for i in range(n_gs):
    best_param[i], _ = optimize_He_gss(He, b_or_z, param_bounds,
                                      energy_args = energy_params,
                                      display = True)
    print('GS run %d completed' % (i+1))
toc = process_time()
print('\nGolden section search time : %5.4f s\n' % (toc - tic))

best_param_med = np.median(best_param) #median is less affected then mean by outliers 

print('Best variational parameter is : %10.7f\n' % best_param_med) 

if b_or_z == 'b' :
    He = Helium(Z=2., b=best_param_med, 
               interacting=interacting, jastrow=jastrow,
               numeric = numeric)
elif b_or_z == 'z':
    He = Helium(Z=best_param_med, b=0., 
               interacting=interacting, jastrow=jastrow,
               numeric = numeric)

print('Evaluating energy with best parameter ...\n')
best_en = np.empty(10)
best_en_err = np.empty(10)
ac_ratio = np.empty(10)
for i in range(9):
    best_en[i], best_en_err[i], ac_ratio[i] = evaluate_energy(He.weight, He.EL, **energy_params)
best_en[-1],  best_en_err[-1], ac_ratio[-1], r6 = evaluate_energy(He.weight, He.EL, **energy_params,
                                                                  save_chain = True, talky = True)

best_en = np.mean(best_en)
best_en_err = np.mean(best_en_err)
ac_ratio = np.mean(ac_ratio)

print('Best value for the Energy functional is ')
print('E = %5.4f +- %5.4f eV$' % (best_en, best_en_err))
print('Acceptance Ratio : %5.4f' % ac_ratio)

print("Saving results into 'Output/He_VMC_results.txt' ")
fh = open('Output/He_VMC_results.txt','w')
fh.write('Best param : %10.8f\n\n' % best_param_med)
fh.write('Best Energy functional value : %10.8f\n' % best_en)
fh.write('Error on the energy functional : %10.8f\n\n' % best_en_err)
fh.write('Metropolis Acceptance ratio on last run : %5.4f' % ac_ratio)
fh.close()

energies = He.EL(r6)
print("Saving chain steps into 'Output/chain.txt'\n")
out = open('Output/chain.txt', 'w')
out.write('     tau         EL         x1         y1         z1         x2         y2         z2\n')
for i in range(0,len(energies[1:]),freq):
    coords = r6[i]
    out.write('%8d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n' 
               % (i, energies[i], coords[0], coords[1], coords[2],
                  coords[3], coords[4], coords[5]) )
out.close()


#graphics 

#Direct search plots
fig, ax = plt.subplots(1,2)
plt.subplots_adjust(wspace=0.5)

_ = ax[0].set_title('Energy functional')
_ = ax[0].set_xlabel(r'Variational parameter $p$')
_ = ax[0].set_ylabel(r'$E[\psi]_p$')

_ = ax[1].set_title('Statistical Error')
_ = ax[1].set_xlabel('Variational parameter $p$')
_ = ax[1].set_ylabel(r'$\Delta E$')

params = np.linspace(param_min, param_max, param_values)
for i in range(n_ds): 
    _ = ax[0].plot(params, ds_energies[i,:], marker='.',
                       alpha=0.5, color='cadetblue', lw=0.5)
    _ = ax[1].plot(params, ds_errors[i,:], marker='.', 
                       alpha=0.5, color='cadetblue', lw=0.5)
    
ax[0].plot(params, avg_ds_energies, c='r', marker='.', label='Average values')
ax[1].plot(params, avg_ds_errors, c='r',marker='.', label='Average values')
fig.savefig('Output/DS_results.png')
print("Saving DS plot in 'Outpu/DS_results.png'")

# Golden section search plots 
fig2, ax2 = plt.subplots()

_ = ax2.set_title('Golden Section Search Results')
_ = ax2.set_xlabel('Gss run index')
_ = ax2.set_ylabel(r'Variational parameter $p$')

_ = ax2.plot(range(n_gs), best_param, color='cadetblue', 
             marker='.', linestyle='none')
_ = ax2.axhline(y=np.median(best_param), label='param median : %5.4f' % np.median(best_param),
                color='r', lw=0.5)
_ = ax2.legend()
fig2.savefig('Output/GSS_results.png')
print("Saving GSS plot in 'Outpu/GSS_results.png'")

#Local Energy histogram
fig3, ax3 = plt.subplots()
_ = ax3.set_title('Local Energy histogram')
_ = ax3.set_xlabel('EL')
_ = ax3.set_ylabel(r'$P(E_L)$')
_ = ax3.hist(energies, bins=200, density=True)
_ = ax3.axvline(x = np.mean(energies), color='r', label = 'Mean energy : %5.4f' % np.mean(energies))
_ = ax3.legend()
fig3.savefig('Output/EL_hist.png')
print("Saving EL histogram in 'Outpu/EL_hist.png'")