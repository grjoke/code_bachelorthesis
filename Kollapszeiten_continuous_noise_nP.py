"""
Python-file for conducting a statistical experiment about the collapse rate of a nonlinear world earth system,
disturbed by continuous noise applied on one variable (P).
It can be called with "python3 Kollapszeiten_continuous_noise_nP.py" .
It generates as well figures as text-file with the results.
"""
import numpy as np
import sdeint as sde #used solver for solving stochastic differential equations
from scipy.integrate import odeint
import sys, os, shutil #librarys for creating, deleting or copying files and directorys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import random as rn
import seaborn as sns
import statsmodels


import parameter as pm #file with all needed parameters for the modell
import Funktion as fc #module for integration of the modell without noise
import Funktion_noiseL_dict as fcnl #module for integration with continuous noise on L
import Funktion_noiseP_dict as fcnp #module for integration with continuous noise on P
import Funktion_LeviNoise as fcln #module for integration with discrete noise, both variables have different functions in this module

#sns.set_theme(color_codes=True)
np.set_printoptions(threshold=sys.maxsize)


#initial values for the simulation:
T = 20000 #time the model should be runned
dt = 1. #stepsize for the integration

#set of initial values, initial value lays in the center of the phasespace
x01 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2880,  # terrestrial carbon [GtC]
 'P': 5e5,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

#set of initial values, initial values lays on the limit cycle of the system
L_init = 2519.5
P_init = 2.64e9
x02 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': L_init,  # terrestrial carbon [GtC]
 'P': P_init,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}


n = 1000 # number of samples

noise_init = 0. #init-value for the noise amplitude
noise_end = 0.1 #end-value for long noise interval

# creating a directory to save all resulting files
if os.path.isdir('files_collaps_rate_gaussian_nP'):
    shutil.rmtree('files_collaps_rate_gaussian_nP')
os.mkdir('files_collaps_rate_gaussian_nP')
#function to save arrays in a textfile in the directory which was created before
def write_results_in_file(array, fname):
    if os.path.exists('files_collaps_rate_gaussian_nP/' + fname + '.txt'):
        os.remove('files_collaps_rate_gaussian_nP/' + fname + '.txt')
    print(array, file = open('files_collaps_rate_gaussian_nP/' + fname + '.txt', 'a'))
    return 0

#function for running the modell with an noise amplitude and calculate the collapsrate, return both arrays at the end
def run_model_with_different_amplitudes(x_init, noise_end):
    # creating the arrays, size = number of samples
    collaps_rate = np.zeros(n)
    amplitude_array = np.zeros(n)
    #run the modell n times, calculate a new noise amplitude in every run
    for i in range(n):
        pm.amplitude_fc = rn.uniform(noise_init, noise_end)
        amplitude_array[i] = pm.amplitude_fc
        sim = fcnp.run_noise(x_init, T, dt, pm.amplitude_fc)
        L = sim['L']/pm.Cges
        for k in range(T):
            if L[k] < 0.01:
                collaps_rate[i] = 1/k
                break
            else:
                continue
    return collaps_rate, amplitude_array

#run the function above for two different sets of initial values and write all return arrays into a file
collaps_rate_1, amplitudes_1 = run_model_with_different_amplitudes(x01, noise_end)
write_results_in_file(collaps_rate_1, fname = 'collaps_rate_random_noise_1_')
write_results_in_file(amplitudes_1, fname = 'amplitudes_collaps_rate_random_noise_1')

collaps_rate_2, amplitudes_2 = run_model_with_different_amplitudes(x02, noise_end)
write_results_in_file(collaps_rate_2, fname = 'collaps_rate_random_noise_2')
write_results_in_file(amplitudes_2, fname = 'amplitudes_collaps_rate_random_noise_2')


#define a linear function for plotting later
def linfit(x, a, b):
    return a*x + b

#calculate a linear regression over the calculated arrays and write the parameters in a file
res1 = linregress(amplitudes_1, collaps_rate_1)
res2 = linregress(amplitudes_2, collaps_rate_2)
write_results_in_file(res1, fname = 'regression_param_collaps_rate_random_1')
write_results_in_file(res2, fname = 'regression_param_collaps_rate_random_2')


#array for plotting smooth functions
amplitude_plot_array = np.linspace(noise_init, noise_end, 1000)
#function for plotting the results with matplotlib
def plot_the_results(amplitude_array, collaps_rate_array, amplitude_plot_array, res_regression, fname):
    sns.regplot(x = amplitude_array, y = collaps_rate_array, robust = False,  label = 'collapse rate', color = 'g')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b', alpha = 1., linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    #plt.plot(label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.ylim(0., 0.002)
    plt.xlabel('continuous noise amplitude')
    plt.ylabel(r'collapse rate $1/T_{collapse}$')
    plt.title('Collapse rate dependent on the noise amplitude \n (continuous noise on P)')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_collaps_rate_gaussian_nP/' + fname + '.pdf')
    plt.savefig('files_collaps_rate_gaussian_nP/' + fname + '.png')
    #plt.show()
    plt.clf()

#running the plot function
plot_the_results(amplitudes_1, collaps_rate_1, amplitude_plot_array, res1, fname = 'collaps_rate_gaussian_noise_1_nP')
plot_the_results(amplitudes_2, collaps_rate_2, amplitude_plot_array, res2, fname = 'collaps_rate_gaussian_noise_2_nP')
