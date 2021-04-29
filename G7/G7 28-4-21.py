# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 12:12:41 2021

@author: Group 7
"""
# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
#import MyTicToc as mt
import pandas as pd
#import PyDREAM as pdr

#%%
#==========================================================================
'Assignment 1'

#definition of parameters
beta0 = 0.85 #0.85
a = 18**1.3 * 24 #*10**-3 #m/day ?? (10^1.3 mm/hr)  initial : 10**1.3 * 24.... first value range: 3to18 
bcl = 7.9 #7.9  upper bound 30
bwb = 28 # upper bound 80
Cf = 2
ncl = 0.30 #porosity for the cover layer
nwb = 0.30 #porosity for the waste body

At = 9100 #m2
Ab = 28355 #m2
#the base area of the cover
H = 1.5 + 12 #m Total height (cover + body)
Acb = At + 1.5/H * (Ab-At)  #m2 (the area of the base of the cover)
Vcl = (At + Acb)/2*1.5 #m3 (volume of the cover)
Vwb = (Ab + Acb)/2*12 #m3 (volume of the waste body)

#for the cover:
Scl_min = 100 #m3, min. water content of the cover layer (cl)
Scl_max = 1000 #ncl*Vcl #1000 #ncl*Vcl #m3, max. water content of the cover layer (cl)

#for the evaporation
Sev_min = 100 #0.001*At #m3 #assumed min evap 
Sev_max = 1000 #0.0057*At #1000  #0.0057*9100 #m3 #assumed max evap from the excel 

#for the waste body
Swb_min = 100 #100 #m3, min. water content of the waste body (wb)
Swb_max = 5000 #nwb*Vwb*10**9 #5000 #nwb*Vwb #m3, maximum water content of the wb (porosity x volume)

#Scl and Swb initial
Scl_0 = 400 #400 #m3/mm3??
Swb_0 = 4000 #m3/mm3??
Y0= np.array([Scl_0, Swb_0])

print('')
print('volume of the cover layer (Vcl) is :', Vcl, 'm3')
print('volume of the waste body layer (Vwb) is :', Vwb, 'm3')
print('')
print('Scl_min is: ', Scl_min, 'm3')
print('Scl_max is: ', Scl_max, 'm3')
print('Sev_min is: ', Sev_min, 'm3')
print('Sev_max is: ', Sev_max, 'm3')
print('Swb_min is: ', Swb_min, 'm3')
print('Swb_max is: ', Swb_max, 'm3')
print('Scl initial is:', Scl_0, 'm3')
print('Swb initial is:', Swb_0, 'm3')
print('')
#%%

# Defining the data from the excel
data = pd.read_excel('G7.xlsx')

'Defining the functions of J, L, and E'
Qdr = np.array(data['Q']) #*10**9 #m3/day or mm3/day
Jrf = np.array(data['Rain']) *1000 #m/day or mm/day     data['Rain'].iloc[:]
pEV = np.array(data['pEV']) *1000 #m/day or mm/day
Temp = data['Temp'].iloc[:] #celcius

#time:
t = np.linspace(0,len(Jrf),len(Jrf)+1) #day

#%%
'Checking the shape of given data'
print('shape of Qdr is:', np.shape(Qdr))
print('shape of Jrf is:', np.shape(Jrf))


#%%
#Define the equation

def dSdt(t, Y, Jrf, pEV):    #Y[0] = Scl, Y[1] = Swb......Y = [Scl, Swb]
    
    #determining Scl, Swb, Lcl, Lwb
    Scl = Y[0][0]
    Swb = Y[1][0]
    Lcl = a*( (Scl - Scl_min)/(Scl_max - Scl_min) )**bcl  #m3
    Lwb = a*( (Swb - Swb_min)/(Swb_max - Swb_min))**bwb  #m3
    
    #determining E & beta
    if Scl < Sev_min:
        fred = 0
    elif Sev_min <= Scl <= Sev_max:
        fred = (Scl - Sev_min)/(Sev_max - Sev_min)
    else: #Scl > Sev_max:
        fred = 1
    
    E = pEV * Cf * fred
    beta = beta0*((Scl - Scl_min)/(Scl_max - Scl_min))
    
    #determining the main equation
    dScl_dt = Jrf - Lcl - E
    dSwb_dt = (1 - beta)*Lcl - Lwb
    
    return dScl_dt, dSwb_dt 

#calculating Qfrom equation 3 
Qdr_new = 0
Qdr_new = [0] #m3 , initial given Q just to make it the same shape as t 
Scl_new = [Scl_0] #an array of Scl_new
Swb_new = [Swb_0] #an array of Swb_new
Q = [beta0 * Scl_0 - Swb_0]

#make a loop because every timestep value J and pEv change which lead to a new integration
for i in range(len(t)-1):
    t_span = [t[i], t[i+1]]
    #solving the ODE equations
    YODE = sp.integrate.solve_ivp(dSdt, t_span, Y0, 
                                   method='RK45', vectorized=True, 
                                   rtol=1e-5, args=(Jrf[i], pEV[i]))
    Y0 = np.array([YODE.y[0][-1], YODE.y[1][-1]])
    Scl_new.append(Y0[0])
    Swb_new.append(Y0[1])
    Qdr_new.append(Qdr[i])
    
#might be handy at a later stage
Scl_new = np.array(Scl_new)
Swb_new = np.array(Swb_new)

#new Lcl and Lwb for next step?
Lcl = a*((Scl_new - Scl_min)/(Scl_max - Scl_min))**bcl
Lwb = a*((Swb_new - Swb_min)/(Swb_max - Swb_min))**bwb
beta = beta0*((Scl_new - Scl_min)/(Scl_max - Scl_min))
Q = (beta*Lcl + Lwb)/1000

Qcum=np.cumsum(Q) #cumulative Outflow
Qcumex=np.cumsum(Qdr_new)

#%%

'Plotting Scl and Swb in the same graph with differrent scale'
#plot the Scl
fig,ax = plt.subplots()
ax.plot(t, Scl_new, color='blue', label='Scl')
ax.set_xlabel('time (d)', fontsize=14)
ax.set_ylabel('Scl (m3)', color='blue',fontsize=14)
plt.legend(bbox_to_anchor=(1.1, 1.0), loc='upper left')
plt.title('Scl and Swb over time',fontsize=15)

#plot the Swb
ax2 = ax.twinx()
ax2.plot(t, Swb_new, color = 'red', label='Swb')
ax2.set_ylabel('Swb (m3)',color='red', fontsize=14)
plt.grid()
plt.legend(bbox_to_anchor=(1.1, 0.9), loc='upper left')
plt.show()

#%%
'Plot Scl vs Swb'
plt.figure()
plt.plot(Scl_new, Swb_new, color='blue', label='Scl vs Swb')
plt.xlabel('Scl_new')
plt.ylabel('Swb_new')
#plt.yscale('log')
plt.title('Scl vs Swb',fontsize=17)
plt.grid()
plt.legend()

#%%
'Plot Q'
plt.figure()
# plt.plot(t, Q, color='black', label='calculated Q')
# plt.plot(t, Qdr_new, color='red', label='given Q')
plt.plot(t, Qcum, color='black', label='calculated Q')
plt.plot(t, Qcumex, color='red', label='given Q')
plt.xlabel('time (d)')
plt.ylabel('Q [m/d]')
plt.title('Q over time (m/d)')
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.grid()

