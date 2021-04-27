#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 09:40:45 2021

@author: jlieftinck
"""


# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import MyTicToc as mt
import pandas as pd


# Definition of parameters

Meteo = pd.read_excel (r'WieringermeerData_Meteo.xlsx')
Leachate = pd.read_excel (r'WieringermeerData_LeachateProduction.xlsx')

rain=Meteo['rain_station'].values
pEV=Meteo['pEV'].values
Temp=Meteo['temp'].values
Meteo_time=Meteo['datetime'].values

LeachateOutflow=Leachate['Outflow'].values # Name column B in Excel 'OutFlow' 
LeachateTime=Leachate['Time'].values


# Calculation of geometries and storage boundarys

Ab = 28355                                      # base area [m²]
At = 9100                                       # top area [m²]
ws = 38                                         # slope width [m]
hcl = 1.5                                       # height cover layer [m]
hwb = 12                                        # height waste body [m]
wwb = 281083000                                 # wet weight waste body [kg]
rhos = 2.65                                     # particle density [t/m³]

Atcl = At                                       # Area on top of the cover layer [m²]
wscl = ws/(hcl+hwb)*hcl                         # slope width along the cover layer[m]
btcl = np.sqrt(Atcl)                            # breadth of cover layer at the top [m]
bbcl = btcl + wscl                              # breadth of cover layer at the bottom [m]
Abcl = bbcl**2                                  # Area of the cover layer at the bottom [m]
Vcl = Atcl * hcl + (Abcl-Atcl) * hcl/2          # Volume of the cover layer [m³]

Atwb = Abcl                                     # Area of the waste body at the top [m²]
Abwb = Ab                                       # Area of the waste body at the bottom [m²]
Vwb = Atwb * hwb + (Abwb-Atwb) * hwb/2          # Volume of the waste body [m³]

rhowb = wwb/1000/Vwb                            # density of waste body [t/m³]
gammawb = wwb/100/Vwb                           # unit weight of waste body [kN/m³]
nwb = (rhos-rhowb)/(rhos-1)                     # porosity of the waste body
Vvwb = nwb*Vwb                                  # total volume of the voids in the body [m³]

rhocl = 1.6                                     # Assumption density cover layer [t/m³]
ncl = (rhos-rhocl)/(rhos-1)                     # porosity cover layer
Vvcl = ncl * Vcl                                # total volume of the voids in the cover layer [m³]

f_red=1
c_f=0.006
betha0 = 1.0                                    # coefficient

Sclmin = 0                                      # minimum storage in cover layer [m³]
Sclmax = Vvcl                                   # maximal storage in cover layer [m³]

Swbmin = 0          # minimum storage in waste body [m³]
Swbmax = Vvwb       # maximum storage in waste body [m³]

a = 10**-6          # hydraulic conductivity [m/d]
bcl = 0.5           # empirical parameter 
bwb = 0.5           # empirical parameter
tm=rain.shape
tm=tm[0]

J=rain
# E=pEV*c_f*f_red

# Lcl = a * ((S[0] - Sclmin)/(Sclmax - Sclmin))**bcl
# Lwb = a * ((S[1] - Swbmin)/(Swbmax - Swbmin))**bwb
# betha = betha0 * ((S[0] - Sclmin) / (Sclmax - Sclmin))
# Qdr = betha * Lcl + Lwb

def g(J,t):
    t=int(t)
    return J[t]

def f(PEv,t):
    t=int(t)
    return PEv[t]

    
def dSdt(t, S):
    """ Return the rate of change for the storage for each layer. """
    return np.array([g - a * ((S[0] - Sclmin) / (Sclmax - Sclmin))**bcl - f*c_f*f_red,
                  (1 - (betha0 * ((S[0]-Sclmin)/(Sclmax - Sclmin)))) * (a * ((S[0] - Sclmin) / (Sclmax - Sclmin))**bcl) - (a * ((S[1] - Swbmin)/(Swbmax - Swbmin))**bwb)])


# def dSdt(t, S):
#     """ Return the rate of change for the storage for each layer. """
#     return np.array([J - a * ((S[0] - Sclmin) / (Sclmax - Sclmin))**bcl - E,
#                   (1 - (betha0 * ((S[0]-Sclmin)/(Sclmax - Sclmin)))) * (a * ((S[0] - Sclmin) / (Sclmax - Sclmin))**bcl) - (a * ((S[1] - Swbmin)/(Swbmax - Swbmin))**bwb)])


# Definition of output times
tOut = np.linspace(0, 6209, 6210)              #time

# Initial case, see table
S0 = np.array([0.01, 0.01])
t_span = [tOut[0], tOut[-1]]
SODE = sp.solve_ivp(dSdt, t_span, S0, t_eval=tOut, 
                          method='RK45', vectorized=True, 
                          rtol=1e-5 )
# infodict['message']                     # >>> 'Integration successful.'
SclODE = SODE.y[0,:]
SwbODE = SODE.y[1,:]

#Qdr=(betha0 * ((S[0]-Sclmin)/(Sclmax - Sclmin)))*(a * ((S[0] - Sclmin) / (Sclmax - Sclmin))**bcl)+(a * ((S[1] - Swbmin)/(Swbmax - Swbmin))**bwb)

# Plot results

# Plot Scl and Swb
plt.figure()
plt.plot(tOut, SclODE, 'r-', label='SclODE')
plt.plot(tOut, SwbODE, 'b-', label='SwbODE')

plt.grid()
plt.legend(loc='best')
plt.xlabel('time')
plt.ylabel('Storage')
plt.title('Evolution')
plt.show()

plt.figure()
plt.plot(SclODE, SwbODE, 'b-', label='ODE')

plt.grid()
plt.legend(loc='best')
plt.xlabel('SwbODE')    
plt.ylabel('SclODE')
plt.title('Evolution')
plt.show()


# Leachate Outflow
plt.figure()
plt.plot(LeachateTime, LeachateOutflow, 'b-', label='Leachate')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Leachate')
plt.show()

# Rainfall
plt.figure()
plt.plot(Meteo_time, rain, 'b-', label='Rainfall')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Rainfall')
plt.show()

# pEV
plt.figure()
plt.plot(Meteo_time, pEV, 'b-', label='pEV')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('pEV')
plt.show()

# Temperature
plt.figure()
plt.plot(Meteo_time, Temp, 'b-', label='Temperature')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Temperature')
plt.show()


