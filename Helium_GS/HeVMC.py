import numpy as np 
import numpy.random as rd 
import numpy.linalg as la 
from scipy.optimize import minimize_scalar
from time import process_time
### Metropolis routines 
#
#
#
def metropolis(steps, init_point, d, weight):
    '''Return a Markov Chain of states evolved through metropolis algorithm.
       
       Arguments:
       -----------------
       steps (integer) : number of final states (except the initial 
       one) returned
       init_point (float,1d array) : initial point in configuration
       space 
       d (float) : semi amplitude of virtual displacement interval
       weight (function) : This function represents the asymptotic
       probability distribution. If properly equilibrated, the chain
       will contain states distributed according to the weight function.
       weight need not to be normalized, but it must work with entries 
       of the same shape as init_point.
       
       Outputs:
       -----------------
       chain ((steps+1 X len(init_point) array): array containing
       the evolved states per rows. 
       acceptance_ratio (float) : ratio equals to the number of accepted moves
       with respect to the total number of moves.'''
    
    dim = len(init_point)
    chain = np.empty((steps+1,dim))
    chain[0] = init_point
    
    #All random extraction are computed in advance 
    nu = rd.randint(low=0, high=dim, size=steps+1)
    disp = rd.uniform(0.,1.,size=steps+1)
    acc_rates = rd.uniform(0.,1.,size=steps+1)
    
    #Evolving chain 
    old_prob = weight(init_point)
    count = 0 
    
    for k in range(1,steps+1):
        # only a coordinate at a time is varied 
        chain[k,:] = chain[k-1,:].copy()
        chain[k,nu[k]] += d * (2.* disp[k] - 1.)
        
        new_prob = weight(chain[k,:]) 
        
        if (new_prob/old_prob) >= acc_rates[k]:
            count += 1 
            old_prob = new_prob 
        else :
            chain[k,nu[k]] = chain[k-1,nu[k]].copy()
        
        #compute acceptance ratio 
        acceptance_ratio = count / steps 
    return chain, acceptance_ratio
#
#
#
# Helium class routines 
def laplacian(f, coords, dx=0.0001):
    '''Compute the Laplacian in cartesian orthogonal coordinates
    in a given point specified by coords.
    As long as the coordinate system is orthogonal, this laplacian
    works in any dimension. 
    '''
    lap = 0.
    
    for i in range(len(coords)):
        c_p = coords.copy()
        c_p[i] += dx 
        c_m = coords.copy()
        c_m[i] -= dx
        lap += ( f(c_p) -2.*f(coords) + f(c_m) ) / dx**2.
    return lap
#
#
#
def LaplacianPsiOverPsi(coords,WaveFunction):
    lap = np.zeros(len(coords))
    delta=0.0001
    for i in range(0,len(coords)):
        tempVal3 = WaveFunction(coords[i,:])
        for j in range(0,len(coords[0])):
            coords[i,j]=coords[i,j]+delta
            tempVal=WaveFunction(coords[i,:])
            coords[i,j]=coords[i,j]-2*delta
            tempVal2=WaveFunction(coords[i,:])
            coords[i,j]=coords[i,j]+delta
            lap[i] += ((tempVal+tempVal2)-2.0*tempVal3 )/ (delta*delta)
        lap[i] /= tempVal3
    return lap
#
#
#
#
#
#
def set_npj_wave_func(Z):
    
    def npj_wave_func(coords): 
        wf = np.exp( -Z*la.norm(coords[:3]) ) * np.exp( -Z*la.norm(coords[3:]) )
        return wf 
    
    return npj_wave_func           
#
#
#
def set_pade_jastrow(b):
    
    def pade_jastrow(coords):
        '''Return the square of the Padé-Jastrow trial wavefunction.
           b is the variational parameter'''
        r1_v = np.array(coords[:3])
        r2_v = np.array(coords[3:])
        r12_v = r1_v - r2_v
        r1 = la.norm(r1_v)
        r2 = la.norm(r2_v)
        r12 = la.norm(r12_v)
        
        
        Z = 2.
        a = 0.5 
        return (np.exp(-Z*r1)*np.exp(-Z*r2)*np.exp((a*r12)/(1.+b*r12))) 
    
    return pade_jastrow
#
#
#
def set_potential_He(interacting):
    
    if interacting==True : 
        
        def potential_He(coords): 
            n = len(coords)
            m = len(coords[0])
            V = np.empty(n)
            Z = 2.
            for i in range(n):
                V[i] = - Z*(1./la.norm( coords[i,:3] ) +\
                            1./la.norm( coords[i,3:] ) ) +\
                            1./la.norm( np.array(coords[i,:3]) - np.array(coords[i,3:]) )
            return V
    
    else :
        
        def potential_He(coords):
            n = len(coords)
            m = len(coords[0])
            V = np.empty(n)
            Z = 2.
            for i in range(n):
                V[i] = - Z*(1./la.norm( coords[i,:3] ) +\
                            1./la.norm( coords[i,3:] ) )
            return V
            
    return potential_He
#
#
#
def set_loc_en_pj(b, interacting):
    
    if interacting == True:
        
        def loc_en_pj(coords):
            '''Return the local energy for the Helium atom, with 
            interacting potential and Padé-Jastrow trial 
            wavefunction. 
            b is the variational parameter'''
            el = np.empty(len(coords))
            for i in range(len(coords)):
                r1_v = np.array(coords[i,:3])
                r2_v = np.array(coords[i,3:])
                r12_v = r1_v - r2_v
                r1 = la.norm(r1_v)
                r2 = la.norm(r2_v)
                r12 = la.norm(r12_v)
    
                Z=2.
                a=0.5
    
                el[i] = -Z**2. + (Z-2.)/r1 + (Z-2.)/r2 +\
                        1./r12 * (1. - 2.*a/((1.+b*r12)**2.) ) +\
                        2.*a*b/((1.+b*r12)**3.) - a**2./((1. + b*r12)**4.) +\
                        Z*a/((1. + b*r12)**2.) * np.dot(r12_v/r12, r1_v/r1 - r2_v/r2)
            return el
    
    else :
        
        def loc_en_pj(coords):
            '''Return the local energy for the Helium atom, with 
            interacting potential and Padé-Jastrow trial 
            wavefunction. 
            b is the variational parameter'''
            el = np.empty(len(coords))
            for i in range(len(coords)):
                r1_v = np.array(coords[i,:3])
                r2_v = np.array(coords[i,3:])
                r12_v = r1_v - r2_v
                r1 = la.norm(r1_v)
                r2 = la.norm(r2_v)
                r12 = la.norm(r12_v)
    
                Z=2.
                a=0.5
    
                el[i] = -Z**2. + (Z-2.)/r1 + (Z-2.)/r2 +\
                        1./r12 * (- 2.*a/((1.+b*r12)**2.) ) +\
                        2.*a*b/((1.+b*r12)**3.) - a**2./((1. + b*r12)**4.) +\
                        Z*a/((1. + b*r12)**2.) * np.dot(r12_v/r12, r1_v/r1 - r2_v/r2) 
            return el
    
    return loc_en_pj
#
#
#
def set_local_energy_npj(Z,interacting=True):
    
    if interacting == True :
        def local_energy(coords):
            el = np.empty(len(coords))
            for i in range(len(coords)):
                el[i] = -Z**2. + (Z-2.)*(1/la.norm(coords[i,:3]) + 1/la.norm(coords[i,3:])) +\
                         1./( la.norm( np.array(coords[i,:3]) - np.array(coords[i,3:]) ) )
            return el
    else :
        def local_energy(coords):
            el = np.empty(len(coords))
            for i in range(len(coords)):
                el[i] = - Z**2. + (Z-2.)*( 1/la.norm(coords[i,:3]) + 1/la.norm(coords[i,3:]) )
            return el
    return local_energy
#
#
# Statistical routines 
def bootstrap_rep(data, func):
    '''Bootstrap replicates'''
    return func(rd.choice(data, len(data)))
#
#
#
def make_stats(chain, pieces, array=False):
    '''Given a chain of data, this is divided into 
    the selected number of subchains. The mean value and the 
    mean value variance are finally returned.'''
    #slice the chain with split
    sub_chains = np.array_split(chain, pieces)
    
    #Evaluate mean values of the sliced pieces 
    means = np.empty(pieces)
    for i in range(pieces):
        means[i] = np.mean(sub_chains[i])
    #return mean value and std deviation
    if array==True:
        return np.mean(means), np.std(means), means 
    else :
        return np.mean(means), np.std(means)
#
#
#
# Energy functional routines 
def evaluate_energy(weight, EL, p0=np.array([1.,1.,1.,-1.,-1.,-1.]),
                    steps=10000, d=1.25, method='bootstrap',
                    samples=1000, talky=False):
    '''This function evaluate the energy functional for a generic system.
    Returns the best value for the functional with a 95% Confidence interval
    around this value.
    Both mean and CI are determined through bootstrap analysis.
    
    Arguments:
    --------------
    weight (func) : weight function for Metropolis
    EL (func) : Local Energy function; must 
    depend only on coordinates. 
    p0 (float, 1d array) : initial point in configuration
    space for metropolis algorithm.
    steps (integer) : number of steps of Markov chain 
    equilibration. Actual run takes 10 times this number 
    of steps.
    d (float) : passed as semi amplitude interval for 
    virtual displacements in metropolis.
    err (bool) : if True error analysis is performed
    bs_samples : number of samples in bootsrap analysis
    talky (bool) : if True error analysis results are printed
    
    Outputs:
    ----------------
    best_energy (float) : Best value for the functional
    CI (float) : confidence interval amplitude; 
    returned only if err == True
    acceptance_ratio (float) : Markov Chain acceptance ratio;
    returned only if err == True
    '''           
    
    r6 = MarkovChain(init_point=p0, prob_dist = weight, d = d) #Markov Chain instance
    r6.set_chain(steps)                                        
    
    #equilibration
    tic = process_time()
    r6.evolve_chain(save=True)
    toc = process_time()
    eq_time = toc - tic
    
    #actual run
    tic = process_time()
    r6.set_chain(10*steps)
    r6.evolve_chain()
    toc = process_time()
    ac_time = toc - tic
    
    #evaluating energy on sampled states
    tic = process_time()
    energies = EL(r6.chain)
    toc = process_time()
    en_time = toc-tic
    
    #evaluating the mìfunctional as a mean value
    tic = process_time()
    if (method.lower() == 'bootstrap') :
        #Bootstrapping
        mean_energies = np.empty(samples)
        #Computing Average energies from bootstrap 
        for i in range(samples):
            mean_energies[i] = bootstrap_rep(energies, np.mean)
        #Taking the mean of the average energies
        best_energy = np.mean(mean_energies)
        percs = np.percentile(mean_energies, [2.5,97.5])
        energy_err = np.abs(percs[0] - percs[1])
     
    elif (method.lower() == 'split') :
        #Splitting MarkovChain in subchains
        best_energy, energy_err, mean_energies = make_stats(energies, samples, array=True)
    else :
        print('Invalid method: you should use bootstrap or split !!!')
        return None
         
    toc = process_time()
    stats_time = toc-tic
        
    if talky==True :
        #detailed error analysis 
        percs = np.percentile(mean_energies, [2.5,97.5])
        std = np.std(mean_energies)
        
        # if you feel talkative this prints a lot of infos ;)
        print('Steps : %d' % (r6.steps))
        print('Acceptance ratio   : %7.6f' % r6.acceptance_ratio)
        print('Average Energy     : %10.9f' % best_energy)
        print('95 % C.I.         : ' , percs)
        print('Standard deviation : %10.9f' % std )
        print('Equilibration time : %7.5f' % eq_time)
        print('Actual run time    : %7.5f' % ac_time)
        print('Energy calc time   : %7.5f' % en_time)
        print('Statistics time     : %7.5f' % stats_time)
        print('************************************************')
            
    return best_energy, energy_err, r6.acceptance_ratio
#
#
#
def optimize_He_ds(Z, b, He, info=False, energy_args=None):
    '''Look for the minimum of the energy functional 
    sampling a finite set of equispaced points in parameter space.
    
    Arguments:
    ------------
    Z (float, 1d array) : Array with Z values 
    b (float, 1d array) : Array with b values
    ! Must have same dimension as Z
    He (Helium object) : Helium object representing the system.
    From this object weight function and local energy are computed.
    err (bool) : if True errors are computed
    
    Output: 
    best_param (float): returns the parameters minimizing the energy.
    best_energy(float): returns the minimum value of the energy.
    errors (float, 1d array) : returned only if err == True.
    Contains the error evaluated at each iteration.
    '''
    energies = np.empty(len(Z))
    errors = np.empty(len(Z))
    
    for i in range(len(Z)):
        He.Z , He.b = Z[i], b[i]
        He.reconfig(He._interacting, He._jastrow)
        
        #check if there are energy_args to pass 
        #into evaluate energy for fine_tuning
        if (type(energy_args).__name__ == 'dict'):
            energies[i], errors[i], _ = evaluate_energy(He.weight, He.EL, He.coords,
                                                        talky=info, **energy_args)
        #if there are not, use default values
        else :
            energies[i], errors[i], _ = evaluate_energy(He.weight, He.EL, He.coords,
                                                        talky=info)
         
    best_indx = np.argmin(energies)
    best_params = (Z[best_indx], b[best_indx])
    best_energy = energies[best_indx]
    
    return best_params, best_energy, energies, errors
#
#
#
def energy_1D(x, He, b_or_z='b', info=False, args=None):
    if b_or_z == 'b':
        He.b = x
    elif b_or_z == 'z':
        He.z = x
    else :
        print('Invalid input:\nyou should put b or z as a string')
        return None
    He.reconfig(He._interacting, He._jastrow)
    
    if (type(args).__name__ == 'dict'):
        energy, _, _ = evaluate_energy(He.weight, He.EL, He.coords,
                                       talky=info, **args)
    else :
        energy, _, _ = evaluate_energy(He.weight, He.EL, He.coords,
                                       talky=info)
    return energy
#
#
#
def optimize_He_gss(He, b_or_z, interval, 
                    info = False, energy_args=None,
                    tol=1e-4, maxit=50, display=False):
    
    #define arguments to pass to energy_1d
    if (type(energy_args).__name__ == 'dict'):
        args = (He, b_or_z, info, energy_args)
    else :
        args = (He, b_or_z, info)

    #define convergence algorithm options
    opts = {'maxiter' : maxit,
            'disp':display,
            'xatol' :tol
           }
    
    #call scipy function minimize_scalar
    #for detailed explanation look at 
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalarhtml#scipy.optimize.minimize_scalar
    optimize_results = minimize_scalar(energy_1D, args=args, tol=None,
                                 bounds=interval, method='bounded',
                                 options=opts)
    best_param = optimize_results.x
    return best_param, optimize_results
#
#
# CLASSES 
class MarkovChain:
    '''Class describing a Markov Chain in phase space,
       starting at init_point and with asymptotic probability 
       distribution given by prob_dist'''
    
    def __init__(self, init_point=np.empty((1,1)), 
                 prob_dist=None, d=0.) :
        self._init_point = init_point     #initial point in phase space 
        self._prob_dist = prob_dist       #asymptotic probability distribution
        self._d = d                       #displacement takes place in [-d,d)
        self.steps = 0                   #Markov chain steps
        self.chain = np.empty(0)         #Store Markov chain entries
        self.acceptance_ratio = 0.       
    
    def set_chain(self, steps):
        if ( len(self._init_point) != 0 ) :
            self.steps = steps
            self.chain = np.empty( (self.steps + 1 , len(self._init_point)) )
            return
        else:
            return
  
    def evolve_chain(self, save=False):
        self.chain, self.acceptance_ratio = metropolis(self.steps,
                                                       self._init_point,
                                                       self._d,
                                                       self._prob_dist)
        if (save == True) : 
            self._init_point = self.chain[-1]
        return
#
#
#
class Helium:
    
    def __init__(self, Z, b, coords=[1.,1.,1.,-1.,-1.,-1.],
                 interacting=True, jastrow=True, numeric=False):
        self.name = 'Helium'
        self.Z = Z                      #Variational parameters Z and b
        self.b = b
        self.coords = coords
        self._interacting = interacting
        self._jastrow = jastrow
        self._numeric = numeric
        self.set_potential()                 
        self.set_WF()
        self.set_weight()
        self.set_local_energy()
    
        
    def set_potential(self):
        self.V = set_potential_He(self._interacting)
        return
    
    
    def set_WF(self):
        if (self._jastrow == True):
            self.WF = set_pade_jastrow(self.b)
        else:
            self.WF = set_npj_wave_func(self.Z)
        return
    
    
    def set_weight(self):
        def weight(coords):
            w = self.WF(coords)**2.
            return w
        self.weight = weight
        return
    
    
    def set_local_energy(self):
        if self._numeric == True :
            #local energy is determined numerically
            def EL(coords):
                el = - 0.5 * LaplacianPsiOverPsi(coords, self.WF)          
                el += self.V(coords)
                return el
        else : 
            #local energy is determined analytically
            if self._jastrow == False:
                EL = set_local_energy_npj(self.Z, self._interacting)
            else:
                EL = set_loc_en_pj(self.b, self._interacting)
        self.EL = EL
        return
    
    def reconfig(self, interacting, jastrow): 
        self._interacting = interacting
        self._jastrow = jastrow
        self.set_potential()                 
        self.set_WF()
        self.set_weight()
        self.set_local_energy()
        return
#
#
#