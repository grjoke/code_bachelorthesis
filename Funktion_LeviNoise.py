import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys
import random as rn

import parameter as pm

np.set_printoptions(threshold=sys.maxsize)


"""
Function for running the system for a certain simulation time with discrete noise applied.
Parameters are the initial condition, the simulation time, the step size and the noise amplitude.
There is one function for both variables L and P, where noise is applied.
These functions use the simulation function of the undisturbed system.
"""

def run(init, tstart, tend, step):
    """
    Run the system from x0 for T years, return values every dt years.

    init is a dictionary with keys 'A', 'S', 'T', ... (containing global values)

    The return value is a similar dictionary,
    its values are arrays whose first axis is the time index,
    it has an additional key 't' for time values.
    """
    def algebra_values(x):
        """
        Extract the initial values and compute the other values from the algebraic
        equations
        """
        # extracting the values from the initial dictionary
        A = x[0]
        L = x[1]
        P = x[2]
        T = x[3]

        # compute the other algebraic values
        M = pm.Cges - L - A #maritime carbon stock
        K = pm.kP*P # whole usable capital

        #XB = pm.a_B*L**2 #auxilliary variable for the energy density of Biomass
        #X = XB #auxilliary variable for all energy densitys
        #factor = (P*K)**0.4/X**0.8 #factor for use of ressources to produce energy
        #B = (XB**0.2 * (P*K)**0.4)/pm.e_B #equation for use of Biomass

        B = pm.b*L**0.4*P**0.6 #alternative equation for biomass use with auxilliary parameter b, doesn't work well,
        E = pm.e_B * B   # total energy output, only biomass is use +

        #Y = pm.y_E * E  #economic production out of biomass energy
        #W = pm.y_E*pm.e_B * B/P + pm.w_L * L/pm.sig #wellbeing because of economic production, which function should i use here??

        #alternative form of the equations for Y and W
        Y = pm.y_B*pm.b*L**0.4*P**0.6
        W = (1-pm.i)*Y/P + (pm.w_L * L)/pm.sig

        #values = [A, B, E, K, L, M, P, T, W, Y]
        #print(values)

        return A, B, E, K, L, M, P, T, W, Y

    def differential_values(logx, unused_t):
        """
        Compute the different flows in the system.
        This function will be used in odeint.
        """
        #get a non-logarithmic array
        x = np.exp(logx)
        #get the values from the algebraic calculations
        A, B, E, K, L, M, P, T, W, Y  = algebra_values(x)  #[v for v in algebra_values(x)]

        #compute the different auxilliary terms
        resp = L*(pm.a_0 + pm.a_T*T)
        photo = L*(pm.l_0 - pm.l_T*T)*np.sqrt(A/pm.sig)
        diff = pm.d * (M - pm.m*A)
        emission = B

        fert = (2*pm.p * W*pm.W_p) / (W**2 + pm.W_p**2)
        mort = (pm.q_0 + pm.q_T*T) / W

        #compute the different flows
        dL_dt = photo - resp - emission #+ rn.gauss(0., 0.5)
        dA_dt = -photo + resp + emission + diff
        dT_dt = pm.g * (A/pm.sig - T) #+ (pm.w_L*L)/pm.sig)   #flow equation for temperature T, without Albedo effect
        dP_dt = P * (fert - mort) #without space competition

        #put the different flows in one array, watch out to have the right order like the initial conditions
        dx_dt = dA_dt, dL_dt, dP_dt, dT_dt

        #print the flow-values in a file, each time the function is called by odeint/sdeint, create a new file and append the values in a new line
        #print(dx_dt/x, file = open("Flows_" + str(pm.a_B) + ".txt", 'a'))

        return dx_dt / x # =d_dx log(x) = 1/x


    # getting the initial conditions and defining the time array
    logx0 = np.log([init['A'], init['L'], init['P'], init['T']])
    t = np.arange(tstart, tend, step)


    #running the whole system with odeint, by giving it the right hand side function, the initial value array an the time-array
    logx, info = odeint(differential_values, logx0, t, full_output = 1)


    # get non logarithmic-array
    x = np.exp(logx)


    #writing the integrated values back in the first function
    A, B, E, K, L, M, P, T, W, Y = algebra_values(x.T)

    #creating the solution dictionary
    return_values = {
        't': t,  # time
        # global state variables:
        'A': A,  # atmospheric carbon [GtC]
        'B': B, # biomass use
        'E': E, # energy used
        'K': K, # sum of capital
        'L': L, # terrestrial carbon [GtC]
        'M': M, # maritime carbon  [GtC]
        'P': P, # Population of Humans [H]
        'T': T,  # global mean surface air temperature [GtC / km²]
        'W': W, # wellbeing
        'Y': Y, # global economic production
    }

    return return_values


def run_LeviNoise_P(init, T_ges, step, amp):
    lambd = 0.001 # variable for exponential-distribution, smaller value --> less impact-points
    keys = ['t', 'A', 'B', 'E', 'K', 'L', 'M', 'P', 'T', 'W', 'Y']
    return_values = {}
    for key in keys:
        return_values.update({key: []})
    T_impact_array = []
    T_lauf = 0
    while T_lauf <= T_ges:
        T_impact = int(np.rint(rn.expovariate(lambd)))
        if T_impact >= T_ges or T_impact == 0:
            continue
        else:
            tstart = int(T_lauf)
            tend = int(T_lauf) + T_impact
            q = rn.uniform(-amp, amp) #variable on how many people die in unexpected event
            sim = run(init, tstart, tend, step) #run the system to the discrete noise event
            A, L, P, T = sim['A'], sim['L'], sim['P'], sim['T']
            for key in keys:
                return_values[key] = np.concatenate((return_values[key], sim[key]))
            #create new init-value dictionary
            init = {
            'A': A[-1],
            'L': L[-1],
            'P': P[-1]*(1-q),
            'T': T[-1]
            }
            T_lauf += T_impact
            T_impact_array.append(T_lauf)
            if T_lauf >= T_ges:
                break
            else:
                continue
    return return_values, T_impact_array

def run_LeviNoise_L(init, T_ges, step, amp):
    lambd = 0.001 # variable for exponential-distribution, smaller value --> less impact-points
    keys = ['t', 'A', 'B', 'E', 'K', 'L', 'M', 'P', 'T', 'W', 'Y']
    return_values = {}
    for key in keys:
        return_values.update({key: []})
    T_impact_array = []
    T_lauf = 0
    while T_lauf <= T_ges:
        T_impact = int(np.rint(rn.expovariate(lambd)))
        if T_impact >= T_ges or T_impact == 0:
            continue
        else:
            #T_impact_array.append(T_impact)
            #if T_lauf + T_impact >= T_ges:
            #    break
            #else:
            tstart = int(T_lauf)
            tend = int(T_lauf) + T_impact
            q = rn.uniform(-amp, amp) #variable on how many forest gets destroyed in unexpected event
            sim = run(init, tstart, tend, step)v #run the system until the discrete noise event
            A, L, P, T = sim['A'], sim['L'], sim['P'], sim['T']
            for key in keys:
                return_values[key] = np.concatenate((return_values[key], sim[key]))
            #create new init-value dictionary
            init = {
            'A': A[-1],
            'L': L[-1]*(1-q),
            'P': P[-1],
            'T': T[-1]
            }
            T_lauf += T_impact
            T_impact_array.append(T_lauf)
            if T_lauf >= T_ges:
                break
            else:
                continue
    return return_values, T_impact_array
