#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 9:25:14 2021

@author: jorislieftinck
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import pandas as pd

#constants
S_cl_min = 100
S_cl_max = 1000
S_wb_min = 100
S_wb_max = 5000
S_Ev_min = 100
S_Ev_max = 1000

a =10**1.3  * 24 #mm/d
b_cl = 7.9
b_wb = 28
beta_0 = 0.85
Cf = 2

#initial conditions
S_cl_start = 400
S_wb_start = 4000
Y0= np.array([S_cl_start,S_wb_start])

Meteo = pd.read_excel (r'WieringermeerData_Meteo.xlsx')
Leachate = pd.read_excel (r'WieringermeerData_LeachateProduction.xlsx')

rain=Meteo['rain_station'].values
J=rain*1000

pEv=Meteo['pEV'].values
pEv=pEv*1000

Temp=Meteo['temp'].values
Meteo_time=Meteo['datetime'].values

#Leachate Excel File
LeachateOutflow=Leachate['Outflow'].values 
LeachateTime=Leachate['Time'].values

#time
t = np.linspace(0,len(J),len(J)+1)

def dYdt(t,Y,J,pEv):
    
    S_cl = Y[0][0]
    S_wb = Y[1][0]
    #print(Y)
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


S_cl_new = [S_cl_start]
S_wb_new = [S_wb_start]
Q = [beta_0*S_cl_start-S_wb_start]

#make a loop because every timestep value J and pEv change which lead to a new integration
for i in range(len(t)-1):
    ts = [t[i],t[i+1]]
    YODE = sp.solve_ivp(dYdt, ts, Y0,
                                   method='RK45', vectorized=True, 
                                   rtol=1e-5, args=(J[i],pEv[i]))
    Y0 = np.array([YODE.y[0][-1],YODE.y[1][-1]])
    S_cl_new.append(Y0[0])
    S_wb_new.append(Y0[1])

#might be handy at a later stage
S_cl_new = np.array(S_cl_new)
S_wb_new = np.array(S_wb_new)

L_cl = a*((S_cl_new-S_cl_min)/(S_cl_max-S_cl_min))**b_cl
L_wb = a*((S_wb_new-S_wb_min)/(S_wb_max-S_wb_min))**b_wb
beta = beta_0*((S_cl_new-S_cl_min)/(S_cl_max-S_cl_min))
Q = (beta*L_cl+L_wb)

Qcum=np.cumsum(Q) #cumulative Outflow

tyear=t*0.00273790926 #days to years

fig,ax = plt.subplots()
ax.plot(tyear,S_cl_new,color='blue', label='shallow water storage')
ax.set_xlabel('time (year)')
ax.set_ylabel('shallow water storage', color='blue',fontsize=14)
plt.title('shallow and deep water storage over time',fontsize=17)
ax2 = ax.twinx()
ax2.plot(tyear,S_wb_new, color = 'red')
ax2.set_ylabel('deep water storage',color='red', fontsize=14)
plt.grid()
plt.show()

plt.plot(tyear,Q,color='black')
plt.xlabel('time (year)')
plt.ylabel('Q [mm/d]')
plt.title('Discharge over time')
plt.grid()

plt.figure()
plt.plot(tyear, Qcum, 'k-', label='Outflow Cumulative')


plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Outflow')
plt.show()


plt.figure()
plt.plot(S_cl_new, S_wb_new, 'y-', label='S_cl vs S_wb')


plt.grid()
plt.legend(loc='best')
plt.xlabel('S_cl')    
plt.ylabel('S_wb')
plt.show()

# Leachate Outflow

Outtime=np.linspace(0,len(LeachateOutflow), len(LeachateOutflow)+1)

Qex=[]
Qexstart=np.array([LeachateOutflow[0]])

for j in range(len(Outtime)-2):
    Qex_j=LeachateOutflow[j+1]-LeachateOutflow[j]
    Qex.append(Qex_j)

Qex=np.array(Qex)
Qex=np.append(Qexstart,Qex)

#LeachateOutflow_check=np.cumsum(Qex)

plt.figure()
plt.plot(LeachateTime, LeachateOutflow, 'b-', label='Leachate Cumulative')


plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Leachate Cumulative')
plt.show()

plt.figure()
plt.plot(LeachateTime, Qex, 'b-', label='Outflow Excel delta')
plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Outflow delta')
plt.show()

# Rainfall
plt.figure()
plt.plot(Meteo_time, rain, 'c-', label='Rainfall')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Rainfall')
plt.show()

# pEV
plt.figure()
plt.plot(Meteo_time, pEv, 'g-', label='pEV')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('pEV')
plt.show()

# Temperature
plt.figure()
plt.plot(Meteo_time, Temp, 'r-', label='Temperature')

plt.grid()
plt.legend(loc='best')
plt.xlabel('Time')    
plt.ylabel('Temperature')
plt.show()

