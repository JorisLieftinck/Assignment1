#!/usr/bin/env python
# coding: utf-8

# In[1]:


# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import MyTicToc as mt
import pandas as pd

# In[40]:


data1 = pd.read_excel('WieringermeerData_LeachateProduction.xlsx')
data2 = pd.read_excel('WieringermeerData_Meteo.xlsx')


# In[41]:


#display(data2)
rain = data2['rain_station'].values
pEV = data2['pEV'].values
temp = data2['temp'].values

# In[42]:

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
beta0 = 1.0                                    # coefficient

Sclmin = 0                                      # minimum storage in cover layer [m³]
Sclmax = Vvcl                                   # maximal storage in cover layer [m³]

Swbmin = 0          # minimum storage in waste body [m³]
Swbmax = Vvwb       # maximum storage in waste body [m³]

a = 10**-6          # hydraulic conductivity [m/d]
bcl = 0.5           # empirical parameter 
bwb = 0.5           # empirical parameter
tm=rain.shape
tm=tm[0]

# In[47]:


# Definition of Rate Equation
def g(rain, t):
    t = int(t)
    return rain[t]
def f(pEV, t):
    t = int(t)
    return pEV[t]

def dSdt(t, Y):
    """ Return the growth rate of fox and rabbit populations. """
    Lcl = a*(Y[0]-Sclmin)/(Sclmax-Sclmin)**bcl

    E= f(pEV, t) * c_f * f_red
    beta = beta0*(Y[0]-Sclmin)/(Sclmax-Sclmin)**bwb
    Lwd = a*(Y[1]-Swbmin)/(Swbmax-Swbmin)

    dscl = g(rain, t) - Lcl-E
    dswb = (1-beta)*Lcl-Lwd
    
    return dscl, dswb


# In[52]:


# Definition of output times
tOut = np.linspace(0, tm-1, tm)               # time

# Initial case, 10 rabbits, 5 foxes
Y0 = np.array([5, 5])
mt.tic()
t_span = [tOut[0], tOut[-1]]
YODE = sp.solve_ivp(dSdt, t_span, Y0, t_eval=tOut, 
                          method='RK45', vectorized=True, 
                          rtol=1e-5)
    
# infodict['message']                     # >>> 'Integration successful.'
rODE = YODE.y[0,:]
fODE = YODE.y[1,:]


# In[53]:


plt.plot(tOut, rODE, 'r-', label='RODE')
plt.plot(tOut, fODE, 'b-', label='FODE')
