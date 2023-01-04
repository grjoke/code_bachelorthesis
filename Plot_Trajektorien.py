"""
Collection of different functions to plot the results of the simulation regarding different variables.
The variables that are plotted are named in the function name.
"""


#importing needed modules and functions from other libraries
import numpy as np
import sdeint as sde
from scipy.integrate import odeint
import sys, os
import matplotlib.pyplot as plt

#importing self written modules out of the same directory
import parameter as pm



def PlotCT(dict, index1, index2, noise):
    plt.semilogy(dict['t'], dict['A'], c = 'g', label = r'$CO_2$ Atmosphäre A')
    plt.semilogy(dict['t'], dict['L'], c = 'c', label = r'$CO_2$ Biomasse L')
    plt.semilogy(dict['t'], dict['M'], c = 'b', linestyle = '-', label = r'$CO_2$ Meere M')
    #plt.semilogy(dict['t'], dict['T'], c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.semilogy(dict['t'], dict['A'] + dict['L'] + dict['M'], c = 'tab:orange', label = r'$CO_2$ gesamt')
    plt.legend(loc = 'best')
    plt.title(r'$CO_2$ und Temperatur')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Kohlenstoff_Temperatur.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Kohlenstoff_Temperatur_noise.pdf')
    plt.clf()
    return 0

def PlotCT_rel(dict, index1, index2, noise):
    plt.plot(dict['t'], dict['A']/pm.Cges, c = 'g', label = r'$CO_2$ Atmosphäre A')
    plt.plot(dict['t'], dict['L']/pm.Cges, c = 'c', label = r'$CO_2$ Biomasse L')
    plt.plot(dict['t'], dict['M']/pm.Cges, c = 'b', linestyle = '-', label = r'$CO_2$ Meere M')
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.plot(dict['t'], (dict['A'] + dict['L'] + dict['M'])/pm.Cges, c = 'tab:orange', label = r'$CO_2$ gesamt')
    plt.legend(loc = 'best')
    plt.title(r'$CO_2$ und Temperatur')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Kohlenstoff_Temperatur_relativ.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Kohlenstoff_Temperatur_relativ_noise.pdf')
    plt.clf()
    return 0

def Plotall(dict, index1, index2, noise):
    plt.semilogy(dict['t'], dict['A'], c = 'g', label = r'$CO_2$ Atmosphäre A')
    #plt.semilogy(dict['t'], dict['B'], c = 'r', label = 'Biomasse B')
    #plt.semilogy(dict['t'], dict['E'], c = 'b', linestyle = '-.', label = 'Energy output E')
    plt.semilogy(dict['t'], dict['K'], c = 'k', label = 'Kapital K')
    plt.semilogy(dict['t'], dict['L'], c = 'c', label = r'$CO_2$ Biomasse L')
    #plt.semilogy(dict['t'], dict['M'], c = 'b', linestyle = '-', label = r'$CO_2$ Meere M')
    plt.semilogy(dict['t'], dict['P'], c = 'y', label = 'Population P')
    plt.semilogy(dict['t'], dict['T'], c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.semilogy(dict['t'], dict['W'], c = 'm', label = 'Wohlstand W')
    plt.semilogy(dict['t'], dict['Y'], c = 'tab:orange', label = 'Economic Production Y')
    plt.legend(loc = 'best')
    plt.title('Alle Trajektorien')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Trajektorien.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Trajektorien_noise.pdf')
    plt.clf()
    return 0

def Plotall_rel(dict, index1, index2, noise):
    plt.plot(dict['t'], dict['A']/pm.Cges, c = 'g', label = r'$CO_2$ Atmosphäre A')
    #plt.semilogy(dict['t'], dict['B'], c = 'r', label = 'Biomasse B')
    #plt.semilogy(dict['t'], dict['E'], c = 'b', linestyle = '-.', label = 'Energy output E')
    #plt.semilogy(dict['t'], dict['K'], c = 'k', label = 'Kapital K')
    plt.plot(dict['t'], dict['L']/pm.Cges, c = 'c', label = r'$CO_2$ Biomasse L')
    plt.plot(dict['t'], dict['M']/pm.Cges, c = 'b', linestyle = '-', label = r'$CO_2$ Meere M')
    plt.plot(dict['t'], dict['P']*1e-9, c = 'y', label = 'Population P')
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.plot(dict['t'], dict['W']/pm.W_p, c = 'm', label = 'Wohlstand W')
    plt.plot(dict['t'], dict['Y']*1e-11, c = 'tab:orange', label = 'Economic Production Y')
    plt.legend(loc = 'best')
    plt.title('Alle Trajektorien')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Trajektorien_relativ.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Trajektorien_relativ_noise.pdf')
    plt.clf()
    return 0

def PlotWP(dict, index1, index2, noise):
    plt.semilogy(dict['t'], dict['K'], c = 'k', label = 'Kapital K')
    plt.semilogy(dict['t'], dict['W'], c = 'm', label = 'Wohlstand W')
    plt.semilogy(dict['t'], dict['Y'], c = 'tab:orange', label = 'Economic Production Y')
    plt.semilogy(dict['t'], dict['P'], c = 'y', label = 'Population P')
    plt.legend(loc = 'best')
    plt.title('Population und Kapital')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Wirtschaft_Population.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Wirtschaft_Population_noise.pdf')
    plt.clf()
    return 0

def PlotWP_rel(dict, index1, index2, noise):
    #plt.plot(dict['t'], dict['K'], c = 'k', label = 'Kapital K')
    plt.plot(dict['t'], dict['W']/pm.W_p, c = 'm', label = 'Wohlstand W')
    plt.plot(dict['t'], dict['Y']*1e-11, c = 'tab:orange', label = 'Economic Production Y')
    plt.plot(dict['t'], dict['P']*1e-9, c = 'y', label = 'Population P')
    plt.legend(loc = 'best')
    plt.title('Population und Kapital')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Wirtschaft_Population_relativ.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Wirtschaft_Population_relativ_noise.pdf')
    plt.clf()
    return 0

def PlotT(dict, index1, index2, noise):
    plt.plot(dict['t'], dict['T'], c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.legend(loc = 'best')
    plt.title('Temperatur T in Kohlenstoffäquivalent')
    plt.xlabel('Zeit in Jahren')
    plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Temperatur.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Temperatur_noise.pdf')
    plt.clf()
    return 0

def Stackplot(dict, index1, index2, noise):
    plt.plot(dict['t'], dict['A']/pm.Cges, c = 'g', label = r'$CO_2$ Atmosphäre A')
    plt.plot(dict['t'], dict['L']/pm.Cges, c = 'c', label = r'$CO_2$ Biomasse L')
    plt.plot(dict['t'], dict['M']/pm.Cges, c = 'b', linestyle = '-', label = r'$CO_2$ Meere M')
    plt.stackplot(dict['t'], dict['A']/pm.Cges, dict['L']/pm.Cges, dict['M']/pm.Cges, baseline = 'zero', colors = ['g', 'c', 'b'])
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = 'Temperatur T')
    plt.plot(dict['t'], (dict['A'] + dict['L'] + dict['M'])/pm.Cges, c = 'tab:orange', label = r'$CO_2$ gesamt')
    plt.legend(loc = 'best')
    plt.title(r'$CO_2$ und Temperatur')
    plt.xlabel('Zeit in Jahren')
    #plt.grid()
    if (noise==False):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Kohlenstoff_Temperatur_relativ_stackplot.pdf')
    elif (noise==True):
        plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Amplitude = '+ str(pm.amplitude_fc) + '/Kohlenstoff_Temperatur_relativ_stackplot_noise.pdf')
    plt.clf()

    return 0

def Plot_Paper(dict, noise):
    plt.plot([], [], c = 'g', label = r'$L/C_{total}$')
    plt.plot([], [], c = 'c', label = r'$A/C_{total}$')
    plt.plot([], [], c = 'b', linestyle = '-', label = r'$M/C_{total}$')
    plt.stackplot(dict['t'], dict['L']/pm.Cges, dict['A']/pm.Cges, dict['M']/pm.Cges, baseline = 'zero', colors = ['g', 'c', 'b'])
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = r'T$\cdot\Sigma/C_{total}$')
    plt.plot(dict['t'], dict['W']/pm.W_p, c = 'm', label = r'W/$W_P$')
    #plt.plot(dict['t'], dict['P']*0.1e-9, c = 'y', label = r'P[$0.1\cdot 10^9$]')
    plt.plot(dict['t'], dict['P']*1e-10, c = 'y', label = r'P[$1\cdot 10^{-10}$]')
    plt.legend(loc = 'best')
    plt.title(r'Overview over all model trajectories')
    plt.xlabel('time in years')
    #plt.axis('off')
    plt.tight_layout()
    #plt.grid()
    if (noise==False):
        #plt.savefig('Plots_wL = ' + str(index1) + '_y_B = ' + str(index2) + '/Trajektorien/Plot_Paper.png')
        plt.savefig('Plot_Paper_short.pdf')
        plt.savefig('Plot_Paper_short.png')
        plt.show()
    elif (noise==True):
        plt.savefig('Plot_Paper_noise.pdf')
        plt.show()
    plt.clf()
    return 0

def Plot_LeviNoise(dict, array):
    plt.plot([], [], c = 'g', label = r'$L/C_{ges}$')
    plt.plot([], [], c = 'c', label = r'$A/C_{ges}$')
    plt.plot([], [], c = 'b', linestyle = '-', label = r'$M/C_{ges}$')
    plt.stackplot(dict['t'], dict['L']/pm.Cges, dict['A']/pm.Cges, dict['M']/pm.Cges, baseline = 'zero', colors = ['g', 'c', 'b'])
    plt.plot(dict['t'], dict['T']/(pm.Cges/pm.sig), c = 'r', linestyle = '--', label = r'T/(C*/$Sigma$)')
    plt.plot(dict['t'], dict['W']/pm.W_p, c = 'm', label = r'W/$W_P$')
    #plt.plot(dict['t'], dict['P']*0.1e-9, c = 'y', label = r'P[$0.1\cdot 10^9$]')
    plt.plot(dict['t'], dict['P']*1e-10, c = 'y', label = r'P[$1\cdot 10^{-10}$]')
    for i in array:
        plt.axvline(x = i, ymin = 0., ymax = 1., color = 'tab:orange', linestyle = ':', alpha = 1.) #label = 'impact times')
    plt.legend(loc = 'best')
    plt.title(r'Overwiev over all model components')
    plt.xlabel('time in years')
    #plt.axis('off')
    plt.tight_layout()
    #plt.grid()
    plt.savefig('Plot_Paper.pdf')
    plt.savefig('Plot_Paper.png')
    plt.show()
    plt.clf()

    return 0
