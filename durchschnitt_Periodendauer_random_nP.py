"""
Python file for conducting a statistic experiment about the average period time of the system with P-noise applied.
It can be called with "python3 durchschnittliche_Periodendauer_random_nP.py" .
It generates as well figures as text-file with the results.
"""


#importing needed modules and functions from other libraries
import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os, shutil
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.special import erf
import random as rn
import seaborn as sns
import statsmodels

#importing self written modules out of the same directory
import parameter as pm
import Funktion as fc
import Funktion_noiseL_dict as fcnl
import Funktion_noiseP_dict as fcnp

#sns.set_theme(color_codes=True)
np.set_printoptions(threshold=sys.maxsize)

#initial values for the simulation
T = 20000 #simulation time
dt = 1. #time step-size

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

n = 1000 # number of samples

noise_init = 0. #lowest value for the noise amplitude
noise_end = 0.1 #maximum value for the noise amplitude

#creating fresh directory to save the results
if os.path.isdir('files_aver_period_time_nP'):
    shutil.rmtree('files_aver_period_time_nP')
os.mkdir('files_aver_period_time_nP')

#function to write the results of a step in a textfile
def write_results_in_file(array, fname):
    if os.path.exists('files_aver_period_time_nP/' + fname + '.txt'):
        os.remove('files_aver_period_time_nP/' + fname + '.txt')
    print(array, file = open('files_aver_period_time_nP/' + fname + '.txt', 'a'))
    return 0

#function for calculating the average period time of the system
def calculate_average_period_time(time_array):
    #if there is only one entry in the time array, then set the values as follows
    if time_array.size <= 1:
        average_period_time = 1.
        number_of_cycles_run = 0.
    #if there are more entries the calculate the time between each entry and then calculate the average mean
    else:
        num = time_array.size
        number_of_cycles_run = num-1
        period_time = np.empty(num-1)
        for i in range(num-1):
            period_time[i] = time_array[i+1] - time_array[i]
        average_period_time = np.sum(period_time)/(num-1)

    return average_period_time, number_of_cycles_run

#function for smoothing the calculated value array and finding the peaks of the population
def find_maxima_in_oscillation(value_array):
    value_array_smooth = savgol_filter(value_array, 1001, 3)
    timepeaks,_ = find_peaks(value_array_smooth, height = 1.28e9, distance = 600) #1000 eventuell zu hoch?
    return timepeaks

#core function of the experiment
def run_model_with_different_amplitudes(x_init):
    #creating arrays to save the results values
    amplitude_array = np.empty(n)
    average_period_time = np.empty(n)
    amplitude_cycles = np.empty(n)
    number_of_cycles = np.empty(n)
    for i in range(n):
        #choosing a noise amplitude and save it
        pm.amplitude_fc = rn.uniform(noise_init, noise_end)
        amplitude_array[i] = pm.amplitude_fc
        amplitude_cycles[i] = pm.amplitude_fc
        #simulate the system and find the peaks, then calculate the average period time
        sim = fcnp.run_noise(x_init, T, dt, pm.amplitude_fc)
        P = sim['P']
        maxP_time = find_maxima_in_oscillation(P)
        average_period_time[i], number_of_cycles[i] = calculate_average_period_time(maxP_time)
        #number_of_cycles[i] = calculate_average_period_time(maxP_time)[1]
    return average_period_time, amplitude_array, number_of_cycles, amplitude_cycles


def clean_period_arrays(value_array, amplitude_array):
    value_array_clean = []
    amplitude_array_clean = []
    n = len(value_array)
    for i in range(n-1):
        if 50 <= value_array[i] <= T:   #< 500
            value_array_clean.append(value_array[i])
            amplitude_array_clean.append(amplitude_array[i])
        else:
            continue
    return value_array_clean, amplitude_array_clean

#run the core function, save the results in different arrays and write these arrays in text files for later use
average_period_time_1, amplitudes_1, number_of_cycles_1, amplitudes_cycles_1 = run_model_with_different_amplitudes(x0)
write_results_in_file(average_period_time_1, fname = 'average_period_time_random_noise_1')
write_results_in_file(amplitudes_1, fname = 'amplitudes_period_time_random_noise_1')
write_results_in_file(number_of_cycles_1, fname = 'number_of_cycles_1')
write_results_in_file(amplitudes_cycles_1, fname = 'amplitudes_cycles_1')

average_period_time_2, amplitudes_2, number_of_cycles_2, amplitudes_cycles_2 = run_model_with_different_amplitudes(x_noise)
write_results_in_file(average_period_time_2, fname = 'average_period_time_random_noise_2')
write_results_in_file(amplitudes_2, fname = 'amplitudes_period_time_random_noise_2')
write_results_in_file(amplitudes_cycles_2, fname = 'amplitudes_cycles_2')

"""
average_period_time_1, amplitudes_1, number_of_cycles_1 = clean_arrays(average_period_time_1, amplitudes_1, number_of_cycles_1)
average_period_time_2, amplitudes_2, number_of_cycles_2 = clean_arrays(average_period_time_2, amplitudes_2, number_of_cycles_2)
"""
average_period_time_1, amplitudes_1 = clean_period_arrays(average_period_time_1, amplitudes_1)
average_period_time_2, amplitudes_2 = clean_period_arrays(average_period_time_2, amplitudes_2)

write_results_in_file(average_period_time_1, fname = 'average_period_time_random_noise_1_clean')
write_results_in_file(amplitudes_1, fname = 'amplitudes_period_time_random_noise_1_clean')

write_results_in_file(average_period_time_2, fname = 'average_period_time_random_noise_2_clean')
write_results_in_file(amplitudes_2, fname = 'amplitudes_period_time_random_noise_2_clean')


#make a linear regression over the average period time and save the calculated parameters in a text file
res1 = linregress(amplitudes_1, average_period_time_1)
res2 = linregress(amplitudes_2, average_period_time_2)
write_results_in_file(res1, fname = 'regression_param_period_time_random_1')
write_results_in_file(res2, fname = 'regression_param_period_time_random_2')

#definition of different functions which can be used for fitting the data for the number of cycles
def linfit(x, a, k):
    return a*x + k

def expfit(x, a, k, cx, cy):
    return a*np.exp(k*x + cx) + cy

def powerlawfit(x, a, b, cy):
    return a*x**b + cy

def hyperbolfit(x, a, k, cx, cy):
    return (a/(x**k-cx)) + cy

def arctanfit(x, a, k, cx, cy):
    return a*np.arctan(k*(x-cx)) + cy

def distfunc_fit(x, a, k, cx, cy):
    return a * erf(k*x + cx) + cy

#function for fitting the results of the measurement of the number of cycles to one of the functions from above
def fit_cycle_array(amplitude_array, value_array):
    popt, pconv = curve_fit(distfunc_fit, amplitude_array, value_array)
    return popt, pconv

#fit the data with a function
popt_1, pconv_1 = fit_cycle_array(amplitudes_cycles_1, number_of_cycles_1)
err1 = np.sqrt(np.diag(pconv_1))
popt_2, pconv_2 = fit_cycle_array(amplitudes_cycles_2, number_of_cycles_2)
err2 = np.sqrt(np.diag(pconv_2))

#write the results of the fitting in a text file
if os.path.exists('files_aver_period_time_nP/param_cyclefit_period_time_random_1.txt'):
    os.remove('files_aver_period_time_nP/param_cyclefit_period_time_random_1.txt')
print('amp = ' + str(popt_1[0]) + ', k = ' + str(popt_1[1]) + ', cx = ' + str(popt_1[2]) + ', cy = ' + str(popt_1[3]) + ', error = ' + str(err1), file = open('files_aver_period_time_nP/param_cyclefit_period_time_random_1.txt', 'a'))
if os.path.exists('files_aver_period_time_nP/param_cyclefit_period_time_random_2.txt'):
    os.remove('files_aver_period_time_nP/param_cyclefit_period_time_random_2.txt')
print('amp = ' + str(popt_2[0]) + ', k = ' + str(popt_2[1]) + ', cx = ' + str(popt_2[2]) + ', cy = ' + str(popt_2[3]) + ', error = ' + str(err2), file = open('files_aver_period_time_nP/param_cyclefit_period_time_random_2.txt', 'a'))
amplitude_plot_array = np.linspace(noise_init, noise_end, 1000)

#function for plotting all the results
def plot_the_results(amplitude_array, period_time_array, amplitude_cycles, cycle_array, res_regression, popt, fname):  #popt,
    plt.scatter(amplitude_array, period_time_array, color = 'r', marker = 'o', label = 'average period time')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b', linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.xlabel('noise amplitude')
    plt.ylabel(r'time in years')
    plt.title('Average period time dependent on the noise amplitude \n (continuous noise on $P$)')
    plt.ylim(500, 3500)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_aver_period_time_nP/' + fname + '.pdf')
    plt.savefig('files_aver_period_time_nP/' + fname + '.png')
    #plt.show()
    plt.clf()

    sns.regplot(x = amplitude_array, y = period_time_array, robust = False, scatter_kws={'alpha':0.3}, label = 'average period time',  marker = 'o', color = 'g')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b',  alpha = 0.7, linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.xlabel('noise amplitude')
    plt.ylabel(r'time in years')
    plt.title('Average period time dependent on the noise amplitude \n (continuous noise on $P$)')
    plt.ylim(500, 3500)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_aver_period_time_nP/' + fname + '_seaborn.pdf')
    plt.savefig('files_aver_period_time_nP/' + fname + '_seaborn.png')
    #plt.show()
    plt.clf()

    plt.scatter(amplitude_cycles, cycle_array, color = 'r', marker = 'o', alpha = 0.5, label = 'number of cycles')
    plt.plot(amplitude_plot_array, distfunc_fit(amplitude_plot_array, *popt), color = 'b', linestyle = '-', label = 'erf-fit')
    plt.xlabel('noise amplitude')
    plt.ylabel('number of cycles')
    plt.title('Number of cycles per run dependent on noise-amplitude \n (continuous noise on $P$)')
    plt.ylim(0.,13)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_aver_period_time_nP/' + fname + '_cycles.pdf')
    plt.savefig('files_aver_period_time_nP/' + fname + '_cycles.png')
    #plt.show()
    plt.clf()



plot_the_results(amplitudes_1, average_period_time_1, amplitudes_cycles_1, number_of_cycles_1, res1, popt_1, fname = 'average_period_time_gaussian_noise_1_nP') #popt_1,
plot_the_results(amplitudes_2, average_period_time_2, amplitudes_cycles_2, number_of_cycles_2, res2, popt_2, fname = 'average_period_time_gaussian_noise_2_nP') #popt_2,
