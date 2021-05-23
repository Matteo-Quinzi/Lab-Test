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
       d (float) : half width of virtual displacement interval
       weight (function) : This function represents the asymptotic
       probability distribution. If properly equilibrated, the chain
       will contain states distributed according to the weight function.
       weight need not to be normalized, but it must work with entries 
       of the same shape as init_point.
       
       Outputs:
       -----------------
       chain ((steps+1 X len(init_point) array): array containing
       the evolved states (one on each row). 
       acceptance_ratio (float) : ratio of the number of accepted moves
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
# INTERNAL HELIUM CLASS ROUTINES
# These are listed out of the Helium class for readability 
# and code maintainability
#
#
def LaplacianPsiOverPsi(coords,WaveFunction):
    ''' Evaluate the Laplacian of a wave function and 
    divide it for the wave function itself.
    
    Arguments: 
    -------------
    coords (2-dim array) : 2-dimensional array containing 
    on each row the coordinates of a state where the wavefunction
    must be evaluated. 
    
    Outputs:
    lap (1-dim array) : 1-dimensional array containing 
    the laplacian of psi over psi evaluated at each state.
    '''
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
    '''Given the variational parameter Z, returns the wavefunction without the 
       Padé-Jastrow approximant.
       The wave function returned works only with 1-dimensional arrays containing
       the two electrons coordinates as [x1,y1,z1,x2,y2,z2]
    '''
    def npj_wave_func(coords): 
        wf = np.exp( -Z*la.norm(coords[:3]) ) * np.exp( -Z*la.norm(coords[3:]) )
        return wf 
    
    return npj_wave_func           
#
#
#
def set_pade_jastrow(b):
    '''Given the variational parameter b, returns the wavefunction without the 
       Padé-Jastrow approximant. Cusp conditions have already been considered. 
       The wave function returned works only with 1-dimensional arrays containing
       the two electrons coordinates as [x1,y1,z1,x2,y2,z2]'''
    def pade_jastrow(coords):
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
    '''Returns the Helium potential. The boolean parameter interacting
       allows to neglect the e-e interaction term.
    '''
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
    '''Returns Helium local energy evaluated for the Padé-Jastrow
       trial wavefunction. b is the variational parameter, while
       the boolean argument interacting allows to neglect e-e interaction 
       term. Cusp conditions are already considered in this implementation.
    '''
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
            non-interacting potential and Padé-Jastrow trial 
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
    '''Returns Helium local energy evaluated for the non-Padé-Jastrow
       trial wavefunction. Z is the variational parameter, while
       the boolean argument interacting allows to neglect e-e interaction 
       term.
    '''
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
# STATISTICAL ROUTINES
def bootstrap_rep(data, func):
    '''Evaluate a bootstrap replicate for a given function.
       
       Arguments:
       --------------
       data (float,array) : array containing the data that must
       be resampled. 
       func (function object) : the replicate is actually the value of the function
       evaluated on the resampled data 
       
       Outputs:
       (float) , returns the replicate evaluated on the resampled data. 
    '''
    return func(rd.choice(data, len(data)))
#
#
#
def make_stats(chain, pieces, array=False):
    '''Given a chain of data, this is divided into 
    the selected number of subchains. The mean value and the 
    mean value standard deviation are finally returned.
    
    Arguments:
    ----------------
    chain (float, array) : the array to be split into subchains;
    pieces (integer) : the number of subchains wanted;
    array (bool) : if True the values of the mean evaluated on each subchain 
    is returned.
    
    Outputs:
    ---------------
    (float)(float) : returns the average of the mean values distribution
    and its standard deviation.
    '''
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
# ENERGY FUNCTIONAL ROUTINES
def evaluate_energy(weight, EL, p0=np.array([1.,1.,1.,-1.,-1.,-1.]),
                    steps=10000, d=1.25, method='bootstrap',
                    samples=1000, talky=False,
                    save_chain=False):
    '''This function evaluate the energy functional for a generic system.
    
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
    samples : number of samples for statistical analysis.
    talky (bool) : if True error analysis results are printed.
    
    Outputs:
    ----------------
    best_energy (float) : Best value for the functional; its evaluated
    as the average of the mean values distribution.
    energy_err (float) : error on the functional value returned;
    if bootstrap analysis is performed it is the 95% C.I. of
    the mean values distribution, otherwise it is the standard deviation.
    acceptance_ratio (float) : Markov Chain acceptance ratio;
    returned only if err == True.
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
            
    if save_chain == False :
        return best_energy, energy_err, r6.acceptance_ratio
    else :
        return best_energy, energy_err, r6.acceptance_ratio, r6.chain
#
#
#
def optimize_He_ds(Z, b, He, info=False, energy_args=None):
    '''Look for the minimum of the energy functional 
    sampling a finite set of points in parameter space.
    
    Arguments:
    ------------
    Z (float, 1d array) : Array with Z values 
    b (float, 1d array) : Array with b values
    ! Must have same dimension as Z
    He (Helium object) : Helium object representing the system.
    From this object weight function and local energy are evaluated.
    info (bool) : if True error analysis is printed out for each run.
    energy_args (dict) : dictionary of key-word arguments to pass to 
    evaluate_energy function. This may be used for fine tuning.
    
    Output: 
    best_param (float,tuple): returns the parameters minimizing the energy.
    best_energy(float): returns the minimum value of the energy.
    energies(float, 1d array) : returns the value of energies evluated on each point.
    errors (float, 1d array) : Contains the error evaluated at each iteration.
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
    '''This function is internally used by optimize_He_gss.
       To use golden section search algorithms a proper
       1-dimensional function must be defined.
       
       Arguments:
       -------------
       x(float) : variable representing the variational parameter.
       He (Helium object) : Helium class object representing the studied
       Helium system.
       b_or_z (character) : determines which variational parameter is being changed;
       info (bool) : if True error analysis results are printed out at each run;
       args (dict) : dictionary containing key-word arguments for evaluate-energy function.
       
       Outputs:
       ---------------
       energy(float) : the best value for the energy functional.
    '''
    if b_or_z == 'b':
        He.b = x
    elif b_or_z == 'z':
        He.Z = x
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
    ''' Function looking for the variational parameter that minimizes 
    Helium energy. Implements a golden section search through scipy 
    function minimize_scalar. Method 'bounded' is chosen as it ensures 
    convergence in the investigated interval and implements 'Brent' algorithm
    to speed up convergence whenever possible. 'Brent' algorithm is based 
    on parabolic interpolation.
    The function must be unstable as the enrgy functional may have considerable 
    fluctuations in the working interval.
    
    Arguments:
    ----------------
    He (Helium object) : Describes the Helium system considered.
    b_or_z (character) : defines which variational parameter is being changed.
    interval (tuple, or list like) : defines the interval boundaries where
    the search is conducted. 'bounded' method ensures convergence within 
    the interval.
    info (bool) : if True error analysis is printed out at each step.
    energy_args (dict) : dictionary of key-word arguments that can be passed to 
    evaluate_energy function for fine-tuning.
    tol (float) : tolerance fr the convergence algorithm. It must be intended
    as absolute tolerance for the 'bounded' algorithm. Note that other methods
    while not implemented ('Brent','golden') use tol as relative tolerance.
    makit(integer) : maximum number of iterations 
    display(bool) : if True a feedback message on convergence success is printed out 
    at the end of the gss algorithm.
    
    Outputs:
    -------------
    best_param (float) : variational parameter that minimizes the functional energy.
    optimize_results (OptimiseResults object) : contains information on the 
    algorithm convergence. More detail on scipy.
    
    '''
    
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
       distribution given by prob_dist.
       
       Attributes: 
       _init_point (float, 1d array) : contains the chain starting 
       point in configuration space.
       _prob_dist (fucn) : gives the asymptotic probability distribution
       of the chain. After proper equilibration, the state of the chains will be 
       distributed with the weight distribution.
       _d (float) : half width of the interval for the virtual displacements in 
       Metropolis algorithm. 
       _steps (integer) : number of states (except the initial one) to sample.
       _chain (2d array): contains the sampled states on each row and 
       the state coordinates on each column.
       acceptance_ratio (float) : number of accepted moves with respet to the 
       total number of moves.
       
       Methods :
       ------------------
       set_chain : initialize an empty Markov chain given the wanted number
       of steps.
       evolve_chain : sample states in configuration space with the weight distribution 
       and using the metropolis function. A boolean parameter 'save' allow to 
       store the last point of the chain as the new init_point (This may be
       helpful when evolving the chain).
    '''
    
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
    ''' Class describing a Helium system. 
    The wavefunction implemented are thought to be used for a 
    variational algorithm minimization procedure.
    In this implementation reduced mass effects are neglected as 
    well as spin effects too.
    This results in having electrons wavefunctions that are symmetric 
    and particles position exchange.
    
    Attributes:
    -----------------
    Z (float) : This is a VARIATIONAL PARAMETER, and not the actual
    charge of Helium that is correctly defined in the potential.
    b (float) : variational parameter used for the Padé-Jastrow trial
    wavefunction.
    coords (float, 1d array) : contains the electrons position with respect
    to the nucleus.
    _interacting(bool) : if True interacting potential is considered.
    _jastrow(bool) : if True Padé-Jastrow trial wavefunction is considered.
    _numeric(bool) : if True the laplacian in the local energy calculation 
    is evaluated numerically. Otherwise the local energy is evaluated through 
    the analytic formula.
    V (func) : evaluate the potential on a 2d array containing the states 
    on each row and the electrons coordinates on each column as 
    [x1,y1,z1,x2,y2,z2].
    EL (func) : evaluate the local energy on a 2d array as described for V.
    WF (func) : evaluate the wavefunction on a 1d array of coordinates.
    weight (func) : evaluate the weight function for Helium as described for WF.
    The weight is easily obtained as the square of WF.
    
    Methods:
    ---------------------
    Note : internal methods preceded by 'set' are just used to 
    set the Helium attributes when a given instance is created.
    
    reconfig : boolean parameters _interacting, _jastrow and _numeric 
    are substitued and the attributes V, EL, WF and weight are 
    reassigned.
    
    '''
    def __init__(self, Z, b, coords=[1.,1.,1.,-1.,-1.,-1.],
                 interacting=True, jastrow=True, numeric=False):
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