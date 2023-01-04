import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os, shutil
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import linregress
import random as rn


import parameter as pm
import Funktion as fc
import Funktion_noiseL_dict as fcnl
import Funktion_noiseP_dict as fcnp

np.set_printoptions(threshold=sys.maxsize)

#initial values for the simulation
T = 20000
dt = 1.

x0 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2880,  # terrestrial carbon [GtC]
 'P': 5e5,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

L_init = 2519.5
P_init = 2.64e9
x_noise = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': L_init,  # terrestrial carbon [GtC]
 'P': P_init,  # human population [humans]
 'T': 800 / pm.sig  # global mean surface air temperature [GtC-equiv.]
}

n = 10 # number of samples

noise_init = 0. #init-value for the noise amplitudegit
noise_end = 0.05

if os.path.isdir('files_period_time'):
    shutil.rmtree('files_period_time')
os.mkdir('files_period_time')

def write_results_in_file(array, fname):
    if os.path.exists('files_period_time/' + fname + '.txt'):
        os.remove('files_period_time/' + fname + '.txt')
    print(array, file = open('files_period_time/' + fname + '.txt', 'a'))
    return 0

"""
def calculate_average_period_time(time_array):
    if time_array.size == 0:
        average_period_time = None
    elif time_array.size == 1:
        num = time_array.size
        period_time = np.empty(num)
        for i in range(num):
            average_period_time = time_array[i]
            #eigentlich nicht werten, auch None zurückgeben
    else:
        num = time_array.size
        period_time = np.empty(num-1)
        for i in range(num-1):
            period_time[i] = time_array[i+1] - time_array[i]
            average_period_time = np.sum(period_time)/(num-1)
    return average_period_time
"""

def calculate_average_period_time(time_array):
    if time_array.size <= 1:
        average_period_time = 1.
        #number_of_cycles_run = 1
    else:
        num = time_array.size
        #number_of_cycles_run = num-1
        period_time = np.empty(num-1)
        for i in range(num-1):
            period_time[i] = time_array[i+1] - time_array[i]
        average_period_time = np.sum(period_time)/(num-1)
    return average_period_time #number_of_cycles_run



def find_maxima_in_oscillation(value_array):
    timepeaks,_ = find_peaks(value_array, height = 1.28e9, distance = 500) #1000 eventuell zu hoch?
    return timepeaks

def run_model_with_different_amplitudes(x_init):
    amplitude_array = np.empty(n)
    average_period_time = np.empty(n)
    for i in range(n):
        pm.amplitude_fc = rn.uniform(noise_init, noise_end)
        amplitude_array[i] = pm.amplitude_fc
        sim = fcnl.run_noise(x_init, T, dt, pm.amplitude_fc)
        P = sim['P']
        maxP_time = find_maxima_in_oscillation(P)
        if maxP_time.size == 0:
            continue
        else:
            average_period_time[i] = calculate_average_period_time(maxP_time)
    return average_period_time, amplitude_array

def clean_value_array(value_array, amplitude_array):
    #return np.array([v for v in value_array if 500 <= v <= 5000])
    value_array_clean = []
    amplitude_array_clean = []
    n = value_array.size
    for i in range(n-1):
        if 500 <= value_array[i] <= 5000:   #< 500
            value_array_clean.append(value_array[i])
            amplitude_array_clean.append(amplitude_array[i])
        else:
            continue
    return value_array_clean, amplitude_array_clean

    """
    for v in value_array:
        if 500 <= v <= 5000:
            value_array_clean.append(v)
            amplitude_array_clean.append(amplitude_array[v])
    return value_array_clean, amplitude_array_clean
    """
    """
    n = value_array.size
    for i in range(n-1):
        if value_array[i] <= 500:   #< 500
            np.delete(value_array, i) #value_array =
            np.delete(amplitude_array, i)
        elif value_array[i] >= 5000:
            np.delete(value_array, i)
            np.delete(amplitude_array, i)
        else:
            continue
    return value_array, amplitude_array
    """


def linfit(x, a, b):
    return a*x + b

average_period_time_1, amplitudes_1 = run_model_with_different_amplitudes(x0)
write_results_in_file(average_period_time_1, fname = 'average_period_time_random_noise_1')
write_results_in_file(amplitudes_1, fname = 'amplitudes_period_time_random_noise_1')

average_period_time_2, amplitudes_2 = run_model_with_different_amplitudes(x_noise)
write_results_in_file(average_period_time_2, fname = 'average_period_time_random_noise_2')
write_results_in_file(amplitudes_2, fname = 'amplitudes_period_time_random_noise_2')


average_period_time_1_clean, amplitudes_1_clean = clean_value_array(average_period_time_1, amplitudes_1)
average_period_time_2_clean, amplitudes_2_clean = clean_value_array(average_period_time_2, amplitudes_2)

write_results_in_file(average_period_time_1_clean, fname = 'period_time_random_noise_1_clean')
write_results_in_file(average_period_time_2_clean, fname = 'period_time_random_noise_2_clean')


res1 = linregress(amplitudes_1_clean, average_period_time_1_clean)
res2 = linregress(amplitudes_2_clean, average_period_time_2_clean)
write_results_in_file(res1, fname = 'regression_param_period_time_random_1')
write_results_in_file(res2, fname = 'regression_param_period_time_random_2')

amplitude_plot_array = np.linspace(noise_init, noise_end, 1000)

def plot_the_results(amplitude_array, period_time_array, res_regression, fname):
    plt.scatter(amplitude_array, period_time_array, color = 'r', marker = 'x', label = 'average period time')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b', linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.xlabel('noise amplitude')
    plt.ylabel(r'average period time in years')
    plt.title('(average) period time dependent on the noise amplitude')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_period_time/' + fname + '.pdf')
    plt.savefig('files_period_time/' + fname + '.png')
    plt.show()
    plt.clf()



plot_the_results(amplitudes_1_clean, average_period_time_1_clean, res1, fname = 'average_period_time_random_noise_1')
plot_the_results(amplitudes_2_clean, average_period_time_2_clean, res2, fname = 'average_period_time_random_noise_2')
