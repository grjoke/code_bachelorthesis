import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os
import matplotlib.pyplot as plt
import random as rn

import parameter as pm
#import Funktion as fc

#parameter values for instant greenhouse-effect and diffussion, leads to other trajectories in the phasespace,
#but stochastic trajectories can't be computed
#pm.d = 1e12
pm.g = 1e12

#define a complex variable
j = 1j
#define the borders of the phaseplot for the different variables
Lmin,Lmax = 1e-5,pm.Cges*(1-1e-5)
Pmin,Pmax = 1,5e9
Amin,Amax = 1e-5,pm.Cges*(1-1e-5)
#unused time variable for function call
t = 10

G = 500 #variable to define the numer of steps in the grid


def phasespaceplot_LP():
    j = 1j
    #create the needed grid for x and y
    ys, xs = np.mgrid[Pmin:Pmax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))

    for i in range(G):
        for j in range(G):
            L = xs[i,j]
            P = ys[i,j]
            A = (pm.Cges-L)/(1+pm.m)
            T = A/pm.sig

            z = [A, L, P, T]
            logz = np.log(z)
            dlogz_dt = differential_values(logz, t)
            dz_dt = dlogz_dt*z

            dA_dt = dz_dt[0]
            dL_dt = dz_dt[1]
            dP_dt = dz_dt[2]
            dT_dt = dz_dt[3]


            dxs[i,j] = dL_dt
            dys[i,j] = dP_dt
            abs = np.sqrt(dxs[i,j]**2 + dys[i,j]**2)

    lw = 2*abs/abs.max()

    plt.streamplot(xs, ys, dxs, dys, density = 1.5, color = 'grey', linewidth = lw)
    plt.xlabel(r'L in GtC')
    plt.ylabel('P in bn')
    plt.title('Phasespace between soil carbon and human population')
    plt.tight_layout()
    #plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) +'/Phasenraum/Phasenraum_LP.pdf')
    plt.savefig('Phasenraum_LP.png')
    plt.savefig('Phasenraum_LP.pdf')
    plt.show()
    plt.clf()
    return 0

def phasespaceplot_LP_noise(index1, index2):
    j = 1j
    #create the needed grid for x and y
    ys, xs = np.mgrid[Pmin:Pmax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))


    for i in range(G):
        for j in range(G):
            L = xs[i,j]
            P = ys[i,j]
            A = (pm.Cges-L)/(1+pm.m)
            T = A/pm.sig

            z = [A, L, P, T]
            logz = np.log(z)
            dlogz_dt = differential_values_noise(logz, t)
            dz_dt = dlogz_dt*z

            dA_dt = dz_dt[0]
            dL_dt = dz_dt[1]
            dP_dt = dz_dt[2]
            dT_dt = dz_dt[3]


            dxs[i,j] = dL_dt
            dys[i,j] = dP_dt
            abs = np.sqrt(dxs[i,j]**2 + dys[i,j]**2)

    lw = 2*abs/abs.max()

    plt.streamplot(xs, ys, dxs, dys, density = 1.5, color = 'grey', linewidth = lw)
    plt.xlabel(r'L in GtC')
    plt.ylabel('P in bn')
    plt.title('phasespace between soil carbon and population\n' + ' with noise (on soil carbon)')
    plt.tight_layout()
    plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) +'/Phasenraum/Phasenraum_LP_noise_amp = ' + str(pm.amplitude_ph) + '.pdf')
    #plt.show()
    plt.clf()
    return 0

def phasespaceplot_AL():
    j = 1j
    #create the needed grid for x and y
    ys, xs = np.mgrid[Amin:Amax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))


    for i in range(G):
        for j in range(G):
            L = xs[i,j]
            P = 1.2e9
            A = ys[i,j]
            T = A/pm.sig
            #Mari = pm.Cges - L - A

            Testwert =  A + L# + Mari
            if Testwert > pm.Cges:
                continue
            else:
                z = [A, L, P, T]
                logz = np.log(z)
                dlogz_dt = differential_values(logz, t)
                dz_dt =  dlogz_dt*z

                dA_dt = dz_dt[0]
                dL_dt = dz_dt[1]
                dP_dt = dz_dt[2]
                dT_dt = dz_dt[3]

                dxs[i,j] = dL_dt
                dys[i,j] = dA_dt

    plt.streamplot(xs, ys, dxs, dys, density = 1.5, color = 'grey') #density = 1.5
    plt.plot(1, 0.42*pm.Cges, marker = 's', markersize = 10, color = 'k', label = '"desert state"')
    plt.plot(3030, 425, marker = 'o', markersize = 10,  color = 'k', label = '"forest state"')
    plt.plot(2405, 625, marker = 'o', markersize = 10, color = 'g',  label = 'unstable equillibrium' )
    plt.xlabel(r'L in GtC')
    plt.ylabel(r'A in GtC')
    plt.title('Phase space between terrestrial carbon and \n'+'atmospheric carbon at fixed population')
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    #plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) +'/Phasenraum/Phasenraum_AL.pdf')
    plt.savefig('Phasenraum_AL.pdf')
    plt.savefig('Phasenraum_AL.png')
    plt.show()
    plt.clf()
    return 0

def phasespaceplot_AL_noise(index1, index2):
    j = 1j
    #create the needed grid for x and y
    ys, xs = np.mgrid[Amin:Amax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))


    for i in range(G):
        for j in range(G):
            L = xs[i,j]
            P = 5e5
            A = ys[i,j]
            T = A/pm.sig
            #Mari = pm.Cges - L - A

            Testwert =  A + L# + Mari
            if Testwert > pm.Cges:
                continue
            else:
                z = [A, L, P, T]
                logz = np.log(z)
                dlogz_dt = differential_values_noise(logz, t)
                dz_dt =  dlogz_dt*z

                dA_dt = dz_dt[0]
                dL_dt = dz_dt[1]
                dP_dt = dz_dt[2]
                dT_dt = dz_dt[3]

                dxs[i,j] = dL_dt
                dys[i,j] = dA_dt

    plt.streamplot(xs, ys, dxs, dys, density = 1.5, color = 'grey')
    plt.xlabel(r'L in GtC')
    plt.ylabel(r'A in GtC')
    plt.title('phasespace between soil carbon and \n'+'atmosphere carbon with fixed population\n'+'and noise on soil carbon')
    plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) +'/Phasenraum/Phasenraum_AL_noise_amp = ' + str(pm.amplitude_ph) + '.pdf')
    #plt.show()
    plt.clf()
    return 0


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
    dL_dt = photo - resp - emission
    dA_dt = -photo + resp + emission + diff
    dT_dt = pm.g * (A/pm.sig - T) #+ (pm.w_L*L)/pm.sig)   #flow equation for temperature T, without Albedo effect
    dP_dt = P * (fert - mort) #without space competition

    #put the different flows in one array, watch out to have the right order like the initial conditions
    dx_dt = dA_dt, dL_dt, dP_dt, dT_dt

    return dx_dt / x               #d_dx log(x) = 1/x

def differential_values_noise(logx, unused_t):
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
    dL_dt = photo - resp - emission + pm.amplitude_ph*rn.normalvariate(0., 0.5)
    dA_dt = -photo + resp + emission + diff #+ pm.amplitude_ph*rn.normalvariate(0., 0.5)
    dT_dt = pm.g * (A/pm.sig - T) #+ (pm.w_L*L)/pm.sig)   #flow equation for temperature T, without Albedo effect
    dP_dt = P * (fert - mort) #+ pm.amplitude_ph*rn.normalvariate(0., 0.5)

    #put the different flows in one array, watch out to have the right order like the initial conditions
    dx_dt = dA_dt, dL_dt, dP_dt, dT_dt

    return dx_dt / x               #d_dx log(x) = 1/x


run = phasespaceplot_AL()
