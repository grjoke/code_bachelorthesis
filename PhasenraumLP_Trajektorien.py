import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os
import matplotlib.pyplot as plt

import parameter as pm
import Funktion as fc
import Funktion_noiseL_dict as fcnl
import Funktion_noiseP_dict as fcnp
#import LeviNoise as ln
import Funktion_LeviNoise as fcln #module for integration with discrete noise, both variables have different functions in this module



#parameter values for instant greenhouse-effect and diffussion, leads to other trajectories in the phasespace,
#but stochastic trajectories can't be computed
#pm.d = 1e12
#pm.g = 1e12

#define a complex variable
j = 1j
#define the borders of the phaseplot for the different variables
Lmin,Lmax = 1e-5,pm.Cges*(1-1e-5)  #2000, 3000
Pmin,Pmax = 1,5e9
#Amin,Amax = 1e-5,pm.Cges*(1-1e-5)
#unused time variable for function call
t_unused = 10

G = 500 #variable to define the numer of steps in the grid

T = 10000
dt = 1.

L0 = 2880
P0 = 5e5
x0 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': L0,  # terrestrial carbon [GtC]
 'P': P0,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}


L_lc = 2350#2795
P_lc = 5.09e7  #1.26e9
x_limit_cycle = {
 'A' : 840,
 'L' : L_lc,
 'P' : P_lc,
 'T' : 800/pm.sig
}

x01 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2534,  # terrestrial carbon [GtC]
 'P': 2.612e9,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

x02 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2530,  # terrestrial carbon [GtC]
 'P': 2.4e9,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

#initial conditions at a critical point to study the noise influence

L_init = 2519.5
P_init = 2.64e9
x_noise = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': L_init,  # terrestrial carbon [GtC]
 'P': P_init,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}


keys = ['t', 'A', 'B', 'E', 'K', 'L', 'M', 'P', 'T', 'W', 'Y']
p = 0.03 #amplitude on event-variable for Levi_noise_Trajectories

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

# compute different trajectories for different initial conditions
"""
traj0 = fcn.run_noise(x0, T, dt, 0.003)
L_sim0 = traj0['L']
P_sim0 = traj0['P']

traj1 = fcn.run_noise(x01, T, dt, 0.003)
L_sim1 = traj1['L']
P_sim1 = traj1['P']

traj2 = fcn.run_noise(x02, T, dt, 0.003)
L_sim2 = traj2['L']
P_sim2 = traj2['P']
"""

"""
n = 4 #number of calculated trajectories, always >= 1 for the odeint-calculated trajectory, must be <= 11
noise_init = 0.002 #init-value for the noise amplitude
noise_step = 0.002 #stepsize-value for the addition of noise to the init value

#defining the arrays to save the values coming from the integration
L_sim = np.zeros((n, T))
P_sim = np.zeros((n, T))
amplitudes = np.zeros(n)

#run one integration without noise
traj = fc.run(x_noise, T, dt)
L_sim[0] = traj['L']
P_sim[0] = traj['P']
amplitudes[0] = 0.
#run several integrations with an enlarging noise-amplitude
for i in range(0, n-1):
    pm.amplitude_fc = noise_init + i*noise_step
    traj = fcnp.run_noise(x_noise, T, dt, pm.amplitude_fc)
    L_sim[i+1] = traj['L']
    P_sim[i+1] = traj['P']
    amplitudes[i+1] = pm.amplitude_fc
"""

def plot_phasenraum():
    j = 1j
    #create the needed grid for x and y
    ys, xs = np.mgrid[Pmin:Pmax:G*j, Lmin:Lmax:G*j]
    dxs = np.zeros((G,G))
    dys = np.zeros((G,G))

    # calculate the flows at every point in the grid for L and P
    for i in range(G):
        for j in range(G):
            #intial values out of the grid
            L = xs[i,j]
            P = ys[i,j]
            A = (pm.Cges-L)/(1+pm.m)
            T = A/pm.sig
            # ćalculate the flows and get them non-logarithmic
            z = [A, L, P, T]
            logz = np.log(z)
            dlogz_dt = differential_values(logz, t_unused)
            dz_dt = dlogz_dt*z
            # write the flows back in an array for every variable
            dA_dt = dz_dt[0]
            dL_dt = dz_dt[1]
            dP_dt = dz_dt[2]
            dT_dt = dz_dt[3]


            dxs[i,j] = dL_dt
            dys[i,j] = dP_dt
            abs = np.sqrt(dxs[i,j]**2 + dys[i,j]**2)

    lw = 2*abs/abs.max()
    return xs, ys, dxs, dys, lw

n = 3 #number of calculated trajectories,

def create_Levi_Noise_trajectories_L(n):
    #defining the arrays to save the values coming from the integration
    L_sim = np.empty((n,T))
    P_sim = np.empty((n,T))
    amplitudes = np.zeros(n)
    for i in range(0, n):
        amplitude = 0.05 + 0.05*i
        amplitudes[i] = amplitude
        traj, impacts = fcln.run_LeviNoise_L(x_noise, T, dt, amplitude)
        for k in range(T):
            L_sim[i,k] = traj['L'][k]
            P_sim[i,k] = traj['P'][k]

    return L_sim, P_sim, amplitudes, n

def create_Levi_Noise_trajectories_P(n):
    #defining the arrays to save the values coming from the integration
    L_sim = np.empty((n,T))
    P_sim = np.empty((n,T))
    amplitudes = np.zeros(n)
    for i in range(0, n):
        amplitude = 0.05 + 0.05*i
        amplitudes[i] = amplitude
        traj, impacts = fcln.run_LeviNoise_P(x_noise, T, dt, amplitude)
        for k in range(T):
            L_sim[i,k] = traj['L'][k]
            P_sim[i,k] = traj['P'][k]

    return L_sim, P_sim, amplitudes, n


traj = fc.run(x_noise, T, dt)
L_sim = traj['L']
P_sim = traj['P']

f = 3200

L_plot = L_sim[f:]
P_plot = P_sim[f:]



L_desert = 40
P_desert = 1e8

L_unstable_left = 2050
P_unstable_left = 5e7

L_unstable_right = 3000
P_unstable_right = P_unstable_left

L_fokus = 2550
P_fokus = 1.25e9


xs, ys, dxs, dys, lw = plot_phasenraum()
#L_sim, P_sim, amplitudes, n = create_Levi_Noise_trajectories_L(n)
#L_sim, P_sim, amplitudes, n = create_Levi_Noise_trajectories_P(n)


colors = ['r', 'b', 'g', 'm', 'c', 'y', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:olive'] #'k',

plt.streamplot(xs, ys, dxs, dys, density = 1.5, color = 'grey', linewidth = lw)
plt.plot(L_plot, P_plot, color = 'k', linewidth = 3, label = 'trajectory of undisturbed system')
plt.plot(L_desert, P_desert, color = 'k', marker = 's', markersize = 10, label = 'desert state fixpoint')
plt.plot(L_fokus, P_fokus, color = 'r', marker = 'o', markersize = 10, label = 'focus fixed point')
plt.plot(L_unstable_left, P_unstable_left, color = 'g', marker = 'o', markersize = 10, label = 'unstable saddle point')
plt.plot(L_unstable_right, P_unstable_right, color = 'g', marker = 'o', markersize = 10, label = 'unstable saddle point')
#plt.plot(L_sim1, P_sim1, color = 'b')
#plt.plot(L_sim2, P_sim2, color = 'g')
#for j in range(0,n):
#    plt.plot(L_sim[j], P_sim[j], color = colors[j], label = r'noise = ' + str(np.round(amplitudes[j]*100)) + r'%')
#plt.plot(L_init, P_init, marker = 'o', markersize = 10, color = 'r', label = 'initial condition 2' )
#plt.plot(L0, P0, marker = 'o', markersize = 10, color = 'k', label = 'initial condition 1' )
plt.plot(L_plot[0], P_plot[0], marker = 'o', markersize = 10, color = 'k', label = 'initial condition' )
#plt.axvline(0.01, ymin=0, ymax = 1, linestyle = ':', linewidth = 5, color = 'tab:orange', label = 'threshold for detecting desert state')
plt.xlabel(r'$L$ in GtC')
plt.ylabel(r'$P$ in bn')
#plt.title('Slice trough phase space along \n terrestrial carbon $L$ and human population $P$ \n' +  'with ' + str(n) + ' exemplary trajectories')
#plt.title('Phasenspace between terrestrial carbon L and human population P \n with both initial conditions')
#plt.title('Slice through the phase space along terrestrial carbon $L$ \n and human population $P$ with both initial conditions \n and the threshold for detecting the desert state')
plt.title('Slice through the phase space along \n terrestrial carbon $L$ and human population $P$ \n with axes $A = C_{total}-L/(1 + m)$ and $T = A/\Sigma$')
plt.legend(loc = 'upper left')
plt.tight_layout()
plt.savefig('Phaseportrait_DiscNoise_L.pdf')
plt.savefig('Phaseportrait_DiscNoise_L.png')
plt.tight_layout()
plt.show()
plt.clf()
