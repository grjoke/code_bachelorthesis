import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os, shutil
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.signal import savgol_filter
from scipy.special import erf
import random as rn
import seaborn as sns
import statsmodels

import parameter as pm
import Funktion as fc
import Funktion_noiseL_dict as fcnl
import Funktion_noiseP_dict as fcnp

np.set_printoptions(threshold=sys.maxsize)
#sns.set_theme(color_codes=True)

#initial values for the simulation
T = 20000
dt = 1.

x01 = {
 'A': 840,  # atmospheric carbon [GtC]
 'L': 2880,  # terrestrial carbon [GtC]
 'P': 5e5,  # human population [humans]
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

n = 1000 # number of samples

noise_init = 0. #init-value for the noise amplitude
noise_end = 0.1

if os.path.isdir('files_period_time_seaborn_nP'):
    shutil.rmtree('files_period_time_seaborn_nP')
os.mkdir('files_period_time_seaborn_nP')

def write_results_in_file(array, fname):
    if os.path.exists('files_period_time_seaborn_nP/' + fname + '.txt'):
        os.remove('files_period_time_seaborn_nP/' + fname + '.txt')
    print(array, file = open('files_period_time_seaborn_nP/' + fname + '.txt', 'a'))
    return 0

def find_maxima_in_oscillation(value_array):
    timepeaks,_ = find_peaks(value_array, height = 1.28e9, distance = 200) #1000 eventuell zu hoch?
    return timepeaks

def calculate_average_period_time(time_array):
    if time_array.size <= 1:
        average_period_time = None
        number_of_cycles_run = None
    else:
        num = time_array.size
        number_of_cycles_run = num-1
        period_time = np.empty(num-1)
        for i in range(num-1):
            period_time[i] = time_array[i+1] - time_array[i]
        average_period_time = np.sum(period_time)/(num-1)

    return average_period_time, number_of_cycles_run

def calculate_period_time(time_array, amplitude):
    if time_array.size <= 1:
        period_time = 1.
        number_of_cycles_run = 0.
    else:
        num = time_array.size
        number_of_cycles_run = num-1
        period_time = np.empty(num-1)
        for i in range(num-1):
            period_time[i] = time_array[i+1] - time_array[i]

    return period_time, number_of_cycles_run


def run_model_with_different_amplitudes(x_init):
    amplitudes_cycle = np.empty(n)
    number_of_cycles = np.empty(n)
    amplitude_array = []
    value_array = []
    #number_of_cycles = []
    for i in range(n):
        pm.amplitude_fc = rn.uniform(noise_init, noise_end)
        sim = fcnp.run_noise(x_init, T, dt, pm.amplitude_fc)
        P_raw = sim['P']
        P_smooth = savgol_filter(P_raw, 1001, 3)
        timepeaks,_ = find_peaks(P_smooth, height = 1.28e9, distance = 600)
        if timepeaks.size <= 1:
            amplitudes_cycle[i] = pm.amplitude_fc
            number_of_cycles[i] = 0.
            #continue
            #amplitude_array.append(pm.amplitude_fc)
            #value_array.append(1.)
            #number_of_cycles.append(1.)
        else:
            num = timepeaks.size
            amplitudes_cycle[i] = pm.amplitude_fc
            number_of_cycles[i] = num-1
            for j in range(num-1):
                amplitude_array.append(pm.amplitude_fc)
                value_array.append(timepeaks[j+1] - timepeaks[j])

    return value_array, amplitude_array, number_of_cycles, amplitudes_cycle

def clean_period_arrays(value_array, amplitude_array): # cycle_array):
    value_array_clean = []
    amplitude_array_clean = []
    n = len(value_array)
    for i in range(n-1):
        if 50 <= value_array[i] <= T:   #< 500
            value_array_clean.append(value_array[i])
            amplitude_array_clean.append(amplitude_array[i])
            #cycle_array_clean.append(cycle_array[i])
        else:
            continue
    return value_array_clean, amplitude_array_clean  #, cycle_array_clean

def clean_cycle_array(cycle_array, amplitudes_cycle):
    cycle_array_clean = []
    amplitudes_cycle_clean = []
    n = len(cycle_array)
    for i in range(n-1):
        if 1 <= cycle_array[i] <= 30:
            cycle_array_clean.append(cycle_array[i])
            amplitudes_cycle_clean.append(amplitudes_cycle[i])
        else:
            continue
    return cycle_array_clean, amplitudes_cycle_clean

period_time_1, amplitudes_1, number_of_cycles_1, amplitudes_cycles_1 = run_model_with_different_amplitudes(x01)
write_results_in_file(period_time_1, fname = 'period_time_random_noise_1')
write_results_in_file(amplitudes_1, fname = 'amplitudes_period_time_random_noise_1')
write_results_in_file(number_of_cycles_1, fname = 'number_of_cycles_1')
write_results_in_file(amplitudes_cycles_1, fname = 'amplitudes_cycles_1')

period_time_2, amplitudes_2, number_of_cycles_2, amplitudes_cycles_2 = run_model_with_different_amplitudes(x02)
write_results_in_file(period_time_2, fname = 'period_time_random_noise_2')
write_results_in_file(amplitudes_2, fname = 'amplitudes_period_time_random_noise_2')
write_results_in_file(number_of_cycles_2, fname = 'number_of_cycles_2')
write_results_in_file(amplitudes_cycles_2, fname = 'amplitudes_cycles_2')

period_time_1, amplitudes_1 = clean_period_arrays(period_time_1, amplitudes_1) #, number_of_cycles_1)
period_time_2, amplitudes_2 = clean_period_arrays(period_time_2, amplitudes_2) # , number_of_cycles_2)

#number_of_cycles_1, amplitudes_cycles_1 = clean_cycle_array(number_of_cycles_1, amplitudes_cycles_1)
#number_of_cycles_2, amplitudes_cycles_2 = clean_cycle_array(number_of_cycles_1, amplitudes_cycles_2)

write_results_in_file(period_time_1, fname = 'period_time_random_noise_1_clean')
write_results_in_file(amplitudes_1, fname = 'amplitudes_period_time_random_noise_1_clean')
#write_results_in_file(number_of_cycles_1, fname = 'number_of_cycles_1_clean')

write_results_in_file(period_time_2, fname = 'period_time_random_noise_2_clean')
write_results_in_file(amplitudes_2, fname = 'amplitudes_period_time_random_noise_2_clean')
#write_results_in_file(number_of_cycles_2, fname = 'number_of_cycles_2_clean')


res1 = linregress(amplitudes_1, period_time_1)
res2 = linregress(amplitudes_2, period_time_2)
write_results_in_file(res1, fname = 'regression_param_period_time_random_1')
write_results_in_file(res2, fname = 'regression_param_period_time_random_2')

def linfit(x, a, k):
    return a*x + k

def expfit(x, a, k, cx, cy):
    return a*np.exp(k*x + cx) + cy

def hyperbolfit(x, a, k, cx, cy):
    return (a/(x**k-cx)) + cy

def powerlawfit(x, a, k, cy):
    return a*x**k + cy

def arctanfit(x, a, k, cx, cy):
    return a*np.arctan(k*(x-cx)) + cy

def distfunc_fit(x, a, k, cx, cy):
    return a * erf(k*x + cx) + cy


#     a * (1 - 0.5*(1 +  erf((k-cx)/(x*np.sqrt(2))))) + cy  #k=theta, cx = mu

def fit_cycle_array(amplitude_array, value_array):
    popt, pconv = curve_fit(distfunc_fit, amplitude_array, value_array) #, p0=[10., 1., 0.03, -1.]
    return popt, pconv

popt_1, pconv_1 = fit_cycle_array(amplitudes_cycles_1, number_of_cycles_1)
err1 = np.sqrt(np.diag(pconv_1))
popt_2, pconv_2 = fit_cycle_array(amplitudes_cycles_2, number_of_cycles_2)
err2 = np.sqrt(np.diag(pconv_2))

if os.path.exists('files_period_time_seaborn_nP/param_cyclefit_period_time_random_1.txt'):
    os.remove('files_period_time_seaborn_nP/param_cyclefit_period_time_random_1.txt')
print('amp = ' + str(popt_1[0]) + ', k = ' + str(popt_1[1]) + ', cx = ' + str(popt_1[2]) + ', cy = ' + str(popt_1[3]) + ', error = ' + str(err1), file = open('files_period_time_seaborn_nP/param_cyclefit_period_time_random_1.txt', 'a')) # ', cx = ' + str(popt_1[2]) + #
if os.path.exists('files_period_time_seaborn_nP/param_cyclefit_period_time_random_2.txt'):
    os.remove('files_period_time_seaborn_nP/param_cyclefit_period_time_random_2.txt')
print('amp = ' + str(popt_2[0]) + ', k = ' + str(popt_2[1]) + ', cx = ' + str(popt_2[2]) + ', cy = ' + str(popt_2[3]) + ', error = ' + str(err2), file = open('files_period_time_seaborn_nP/param_cyclefit_period_time_random_2.txt', 'a')) #', cx = ' + str(popt_2[2]) + #', cx = ' + str(popt_2[2]) + ', cy = ' + str(popt_2[3]) +

amplitude_plot_array = np.linspace(noise_init, noise_end, 1000)

def plot_the_results(amplitude_array, period_time_array, cycle_array, amplitudes_cycle, res_regression, popt, fname): #popt
    plt.scatter(amplitude_array, period_time_array, color = 'r', marker = 'o', alpha = 0.3, label = 'period times')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b', linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.xlabel('noise amplitude')
    plt.ylabel(r'time in years')
    plt.title('Cycle duration dependent on the noise amplitude \n (continuous noise on $P$)')
    plt.ylim(500, 5000)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_period_time_seaborn_nP/' + fname + '.pdf')
    plt.savefig('files_period_time_seaborn_nP/'+ fname + '.png')
    #plt.show()
    plt.clf()

    sns.regplot(x = amplitude_array, y = period_time_array, robust = False, scatter_kws={'alpha':0.3}, label = 'period time', marker = 'o', color = 'g')
    plt.plot(amplitude_plot_array, linfit(amplitude_plot_array, res_regression.slope, res_regression.intercept), color = 'b',  alpha = 0.5, linestyle = '-', label = 'linear regression, slope = ' + str(np.round(res_regression.slope, 4)))
    plt.xlabel('noise amplitude')
    plt.ylabel(r'time in years')
    plt.title('Cycle duration dependent on the noise amplitude \n (continuous noise on $P$)')
    plt.ylim(500, 5000)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_period_time_seaborn_nP/' + fname + '_seaborn.pdf')
    plt.savefig('files_period_time_seaborn_nP/'+ fname + '_seaborn.png')
    #plt.show()
    plt.clf()

    plt.scatter(amplitudes_cycle, cycle_array, color = 'r', marker = 'o', alpha = 0.5, label = 'number of cycles')
    plt.plot(amplitude_plot_array, distfunc_fit(amplitude_plot_array, *popt), color = 'b', linestyle = '-', label = 'erf-fit')
    plt.xlabel('noise amplitude')
    plt.ylabel('number of cycles')
    plt.title('Number of cycles per run dependent on noise-amplitude \n (continuous noise on $P$)')
    plt.ylim(0.,13)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('files_period_time_seaborn_nP/' + fname + '_cycles.pdf')
    plt.savefig('files_period_time_seaborn_nP/' + fname + '_cycles.png')
    #plt.show()
    plt.clf()



plot_the_results(amplitudes_1, period_time_1, number_of_cycles_1, amplitudes_cycles_1, res1, popt_1, fname = 'period_time_gaussian_noise_1_nP') # popt_1,
plot_the_results(amplitudes_2, period_time_2, number_of_cycles_2, amplitudes_cycles_2, res2, popt_2, fname = 'period_time_gaussian_noise_2_nP') # popt_2,
