import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys
import random as rn

import parameter as pm

np.set_printoptions(threshold=sys.maxsize)

"""
Function for simulating the undisturbed system.
Parameters are the simulation time and the step size as well as the initial condition.
The function consists of two inner functions where the algebraic as well as the differential equations are implemented.
"""

def run(init, Time, step):
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
        E = pm.e_B * B   # total energy output, only biomass is used

        #Y = pm.y_E * E  #economic production out of biomass energy
        #W = pm.y_E*pm.e_B * B/P + pm.w_L * L/pm.sig #wellbeing because of economic production, which function should i use here??

        #alternative form of the equations for Y and W
        Y = pm.y_B*pm.b*L**0.4*P**0.6
        W = (1-pm.i)*Y/P + (pm.w_L * L)/pm.sig

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
        dL_dt = photo - resp - emission
        dA_dt = -photo + resp + emission + diff
        dT_dt = pm.g * (A/pm.sig - T) #+ (pm.w_L*L)/pm.sig)   #flow equation for temperature T, without Albedo effect
        dP_dt = P * (fert - mort) #without space competition

        #put the different flows in one array, watch out to have the right order like the initial conditions
        dx_dt = dA_dt, dL_dt, dP_dt, dT_dt

        return dx_dt / x # =d_dx log(x) = 1/x


    # getting the initial conditions and defining the time array
    logx0 = np.log([init['A'], init['L'], init['P'], init['T']])
    t = np.arange(0, Time, step)


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
