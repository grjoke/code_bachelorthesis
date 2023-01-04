import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, shutil, os
import matplotlib.pyplot as plt



import parameter as pm
import Funktion as fc
#import Funktion_noise_array as fcn
import Funktion_noiseL_dict as fcnl
import Funktion_noiseP_dict as fcnp

import Plot_Trajektorien as tr
#import Plot_Phasenraum as ph
#import PhasenraumLP_Trajektorien as ph2
#import runandplot as rp

np.set_printoptions(threshold=sys.maxsize)

"""
# parameter values for instant greenhouse-effect and diffussion, leads to other trajectories in the phasespace,
#but stochastic trajectories can't be computed
pm.d = 1e12
pm.g = 1e12
"""

#define important parameters, take the oscillating model
pm.w_L = 4.425e7
pm.y_B = 2.47e9

P0= 5e5 #inital value for human population

T = 10000#Time the system should be runned in years
dt = 1 # Value of stepsize for the timearray


#set initial state for global variables, taken from Nitzborn, Heitzig et. al.,
x0 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2880,  # terrestrial carbon [GtC]
 'P': P0,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

x01 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2534,  # terrestrial carbon [GtC]
 'P': 2.612e9,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

L_init = 2519.5
P_init = 2.64e9
x02 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': L_init,  # terrestrial carbon [GtC]
 'P': P_init,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

x0array = [840, 2880, P0, 800/pm.sig]

pm.amplitude_fc = 0.
pm.amplitude_ph = 1.

run = fcnl.run_noise(x0, T, dt, 0.01)
#print(np.size(run))
#print(run, file = open('Werte.txt', 'a'))
tr.Plot_Paper(run, noise = True)

"""
def Plot_Paper(dict, noise):
    plt.plot([], [], c = 'g', label = r'$L/C_{ges}$')
    plt.plot([], [], c = 'c', label = r'$A/C_{ges}$')
    plt.plot([], [], c = 'b', linestyle = '-', label = r'$M/C_{ges}$')
    plt.stackplot(dict['t'], dict['L']/pm.Cges, dict['A']/pm.Cges, dict['M']/pm.Cges, baseline = 'zero', colors = ['g', 'c', 'b'])
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = r'T/(C*/$Sigma$)')
    plt.plot(dict['t'], dict['W']/pm.W_p, c = 'm', label = r'W/$W_P$')
    #plt.plot(dict['t'], dict['P']*0.1e-9, c = 'y', label = r'P[$0.1\cdot 10^9$]')
    plt.plot(dict['t'], dict['P']*1e-10, c = 'y', label = r'P[$1\cdot 10^{-10}$]')
    plt.legend(loc = 'best')
    plt.title(r'Overwiev over all model components')
    plt.xlabel('time in years')
    #plt.axis('off')
    plt.tight_layout()
    #plt.grid()
    if (noise==False):
        #plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Plot_Paper.png')
        plt.savefig('Plot_Paper.pdf')
        plt.savefig('Plot_Paper.png')
        plt.show()
    elif (noise==True):
        plt.savefig('Plot_Paper_noise.pdf')
        plt.savefig('Plot_Paper_noise.png')
        plt.show()
    plt.clf()

plot = Plot_Paper(run, noise = False)
"""
#run_s = fcn.run_noise(x0array, T, dt, pm.amplitude_fc)

#plot = rp.runplot_without_noise(x0, T, dt, pm.w_L, pm.y_B)
