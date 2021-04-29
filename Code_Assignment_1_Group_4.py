#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 9:25:14 2021

@author: Group 4
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import pandas as pd

# Constants

Ab = 28355                                      # base are [m²]
At = 9100                                       # top area [m²]
ws = 38                                         # slope width [m]
hcl = 1.5                                       # height cover layer [m]
hwb = 12                                        # height waste body [m]
wwb = 281083000                                 # wet weight waste body [kg]
rhos = 2.65                                     # particle density [t/m³]

# Calculating the Storage values
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

S_cl_max = Vvcl
S_wb_max = Vvwb
S_cl_min = 100.0
S_wb_min = 1000.0

#Assumed values for Storage

S_cl_min = 100  #m3
S_cl_max = 1000 #m3
S_wb_min = 100  #m3
S_wb_max = 5000 #m3
S_Ev_min = 100  #m3
S_Ev_max = 1000 #m3

a=18**0.8*24 
b_cl = 7.9
b_wb = 28
beta_0 = 0.8
Cf = 2

# Initial conditions

S_cl_start = 400    #m3
S_wb_start = 4000   #m3
Y0= np.array([S_cl_start,S_wb_start])

# Excel Data (Data has been cut from 2012)

Meteo = pd.read_excel (r'WieringermeerData_Meteo.xlsx')
Leachate = pd.read_excel (r'WieringermeerData_LeachateProduction.xlsx')
rain=Meteo['rain_station'].values
J=rain*1000
pEv=Meteo['pEV'].values
pEv=pEv*1000
Temp=Meteo['temp'].values
Meteo_time=Meteo['datetime'].values

# Cumulative Leachate Excel File

LeachateOutflow=Leachate['Outflow'].values  #Named the 2nd column Outflow in Excel
LeachateOutflow=LeachateOutflow / Ab #in m/day
LeachateTime=Leachate['Time'].values

# Timespan

t = np.linspace(0,len(J),len(J)+1)

# ODE Equation 

def dYdt(t,Y,J,pEv):
    
    S_cl = Y[0][0]
    S_wb = Y[1][0]
    
    L_cl = a*((S_cl-S_cl_min)/(S_cl_max-S_cl_min))**b_cl
    L_wb = a*((S_wb-S_wb_min)/(S_wb_max-S_wb_min))**b_wb
    
    if S_cl < S_Ev_min:
        f_red = 0
    elif S_Ev_min <= S_cl <= S_Ev_max:
        f_red = (S_cl-S_Ev_min)/(S_Ev_max-S_Ev_min)
    elif S_cl > S_Ev_max:
        f_red = 1

    E = pEv*Cf*f_red
    beta = beta_0*((S_cl-S_cl_min)/(S_cl_max-S_cl_min))
    
    #final equations
    dS_cl = J-L_cl-E
    dS_wb = (1-beta)*L_cl-L_wb
    return dS_cl, dS_wb

Qdr_new = [0]  #m3 , initial given Q, same shape as t 
S_cl_new = [S_cl_start]
S_wb_new = [S_wb_start]
Q = [beta_0*S_cl_start-S_wb_start]

#Loop for J and pEv and the ODE solver

for i in range(len(t)-1):
    ts = [t[i],t[i+1]]
    YODE = sp.solve_ivp(dYdt, ts, Y0,
                                   method='RK45', vectorized=True, 
                                   rtol=1e-5, args=(J[i],pEv[i]))
    Y0 = np.array([YODE.y[0][-1],YODE.y[1][-1]])
    S_cl_new.append(Y0[0])
    S_wb_new.append(Y0[1])
    Qdr_new.append(LeachateOutflow[i]) #excel data
    
# Define the values    

S_cl_new = np.array(S_cl_new)
S_wb_new = np.array(S_wb_new)

L_cl = a*((S_cl_new-S_cl_min)/(S_cl_max-S_cl_min))**b_cl
L_wb = a*((S_wb_new-S_wb_min)/(S_wb_max-S_wb_min))**b_wb
beta = beta_0*((S_cl_new-S_cl_min)/(S_cl_max-S_cl_min))
Q = (beta*L_cl+L_wb)/1000

Qcum=np.cumsum(Q) #cumulative Outflow
tyear=t/365 #days to years

# Non-cumulative outflow Excel

Outtime=np.linspace(0,len(Qdr_new), len(Qdr_new)+1)

Qex=[]
Qexstart=np.array([Qdr_new[0]])

for j in range(len(Outtime)-2):
    Qex_j=Qdr_new[j+1]-Qdr_new[j]
    Qex.append(Qex_j)

Qex=np.array(Qex)
Qex=np.append(Qexstart,Qex)


# Storage S_cl and S_wb

fig,ax = plt.subplots()
ax.plot(tyear,S_cl_new,color='blue')
ax.set_xlabel('time (year)')
ax.set_ylabel('Storage Cover Layer [m3]', color='blue',fontsize=14)
plt.title('Storage Cover and Wase Layer',fontsize=17)
ax2 = ax.twinx()
ax2.plot(tyear,S_wb_new, color = 'red')
ax2.set_ylabel('Storage Waste layer [m3]',color='red', fontsize=14)
plt.grid()
plt.show()


# Leachate Outflow Q, cumulative 

plt.figure()
plt.plot(tyear, Qcum, 'k-', label='Calculated Outflow rate Cumulative')
plt.plot(tyear, Qdr_new, 'g-', label='Given Outflow rate Cumulative')
plt.grid()
plt.legend(loc='best')
plt.xlabel('Time [year]')    
plt.ylabel('Outflow Q [m/day]')
plt.title('Discharge cumulative over time')
plt.show()

# Non-cumulative plot Q, excel

plt.figure()
plt.plot(tyear, Qex, 'g-', label='Outflow Excel, delta')

# Non-Cumulative plots Q, calculated

plt.plot(tyear,Q, 'k-', label='Outflow Calculated, delta')
plt.xlabel('Time [year]')
plt.ylabel('Outflow Q [m/day]')
plt.title('Discharge delta over time')

plt.grid()
plt.legend(loc='best')
plt.show()




# # Combi storage plots
 
# plt.figure()
# plt.plot(S_cl_new, S_wb_new, 'y-', label='S_cl vs S_wb')


# plt.grid()
# plt.legend(loc='best')
# plt.xlabel('S_cl')    
# plt.ylabel('S_wb')
# plt.show()

# # Cumulative Storage plots

# Scl=np.cumsum(S_cl_new)
# Swb=np.cumsum(S_wb_new)

# plt.figure()
# plt.plot(t, Scl, 'k-', label='Storage Cumulative')
# plt.plot(t, Swb, 'k-', label='Storage Cumulative')

# plt.grid()
# plt.legend(loc='best')
# plt.xlabel('Time')    
# plt.ylabel('Storage')
# plt.show()

# # Rainfall
# plt.figure()
# plt.plot(Meteo_time, rain, 'c-', label='Rainfall')

# plt.grid()
# plt.legend(loc='best')
# plt.xlabel('Time')    
# plt.ylabel('Rainfall')
# plt.show()

# # pEV
# plt.figure()
# plt.plot(Meteo_time, pEv, 'g-', label='pEV')

# plt.grid()
# plt.legend(loc='best')
# plt.xlabel('Time')    
# plt.ylabel('pEV')
# plt.show()

# # Temperature
# plt.figure()
# plt.plot(Meteo_time, Temp, 'r-', label='Temperature')

# plt.grid()
# plt.legend(loc='best')
# plt.xlabel('Time')    
# plt.ylabel('Temperature')
# plt.show()


