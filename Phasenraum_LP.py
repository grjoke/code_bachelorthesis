import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os
import fileinput
import matplotlib.pyplot as plt

import parameter as pm
#import Funktion as fc

np.set_printoptions(threshold=sys.maxsize)

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
    W = pm.y_B*B/P + (pm.w_L * L)/pm.sig

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
    dL_dt = photo - resp - emission
    dA_dt = -photo + resp + emission + diff
    dT_dt = pm.g * (A/pm.sig - T) #+ (pm.w_L*L)/pm.sig)   #flow equation for temperature T, without Albedo effect
    dP_dt = P * (fert - mort) #without space competition

    #put the different flows in one array, watch out to have the right order like the initial conditions
    dx_dt = dA_dt, dL_dt, dP_dt, dT_dt

    return dx_dt / x              #d_dx log(x) = 1/x


#define a complex variable
j = 1j
#define the borders of the phaseplot and
Lmin,Lmax = 1e-10,1.
Pmin,Pmax = 1e-10,1.5
#unused time variable for function call
t = 10

G = 500 #variable to define the numer of steps in the grid

def phasespaceplot_LP(index):
    #create the needed grid for x and y
    ys, xs = np.mgrid[Pmin:Pmax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))

    for i in range(G):
        for j in range(G):
            L = xs[i,j]*(pm.Cges-1)
            P = ys[i,j]*1e15
            A = (pm.Cges-L)/(1+pm.m)
            T = A/pm.sig

            z = [A, L, P, T]
            logz = np.log(z)
            dlogz_dt = differential_values(logz, t)
            dz_dt = dlogz_dt#*z

            dA_dt = dz_dt[0]
            dL_dt = dz_dt[1]
            dP_dt = dz_dt[2]
            dT_dt = dz_dt[3]

            dxs[i,j] = dL_dt
            dys[i,j] = dP_dt



    plt.streamplot(xs, ys, dxs, dys)
    plt.xlabel(r'L/$C_{ges}$')
    plt.ylabel('P in bn')
    plt.title('phasespace between soil carbon and population')
    plt.savefig('Plots_wL = ' + str(index) + '/Phasenraum_LP.pdf')
    #plt.show()
    plt.clf()
    return 0
