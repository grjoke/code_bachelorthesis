"""
Collection of parameters used in the model.
All used parameters are defined and explained here and are then used in the code by impporting the module with
    import parameter as pm
"""
import numpy as np
import math as m


amplitude_ph = 5. #amplitude for the noise-distribution in creating the phasespaceplots
amplitude_fc = 0.001 #amplitude of the noise distribution in the modell function

#a stands for years, so a^-1 means per year

sig = 1.5e8  #Available Earth Surface in km^2
Cges = 4000 # if Pre-Industrial (no geological Carbon)/5500 if with geologiocal carbon  #total available Carbon Stock in GtC
Cges_PI = 4000 #Pre-Industrial total Carbon stock in GtC
a_0 = 0.0298 #respiration baseline coefficient in a^-1
a_T = 3.2e3 #Respiration sensitivity to temperature in km^2 a^-1 GtC^-1
l_0 = 26.4 #Photosynthesis baseline coefficient in km a^-1 GtC^-1/2
l_T = 1.1e6 #Photosynthesis sensitivity to temperature in km^3 a^-1 GtC^-3/2
d = 0.016 #diffusion rate in a^-1, eventuell sehr groß einstellen, instantane Diffusion
m = 1.43 #solubility coefficient as a number

g = 0.02 #speed of greenhouse effect, factor in diff. equation, eventuell sehr groß einstellen, instantaner Treibhauseffekt


"""
# parameter values for instant greenhouse-effect and diffussion, leads to other trajectories in the phasespace,
#but stochastic trajectories can't be computed
d = 1e12
g = 1e12
"""

p = 0.04 #fertility rate maximum in a^-1
W_p = 2000  #fertility saturation well-being in a^-1 H^-1
q_0 = 20  #mortality baseline coefficient in a^-2 H^-1
q_T = 0   #temperature dependent mortality coefficient
q_P = 0  #mortality coefficient because of space space competition

i = 0 #0.25 #investment ratio
k = 0.1 #capital deprecation rate in a^1

kP = 1000 #fixed capital per Person

a_B = 1  #biomass sector productivity, can be varied
a_F = 0  #fossil fuel sector productivity, can be varied, is zero because we observe a precapitalistic modell
e_B = 40e9 #Biomass energy density in GJ GtC^-1
e_F = 4e10 # fossil fuel energy density in GJ GtC^-1
y_E = 147 #economic output per energy input in $ GJ^-1

b = 5.4e-7 #biomass harvesting rate in GtC^3/5 a^-1 H^-3/5
#b = (a_B**0.2*kP**0.4)/e_B #auxilliary parameter, which combines different parameters for the use of biomass

w_L = 4.425e7 #3.47e7 #well-being sensitivity to land Carbon, can be varied
y_B = 2.47e9 #economic output per biomass input
#y_B = y_E*e_B
