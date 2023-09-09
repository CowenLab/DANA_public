# -*- coding: utf-8 -*-
"""
Montague model.
Montague PR, McClure SM, Baldwin PR, Phillips PEM, Budygin E a, Stuber GD, et al. Dynamic gain control of dopamine delivery in freely moving animals. The Journal of Neuroscience : The Official Journal of the Society for Neuroscience. 2004;24:1754–1759.

"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Use the values from a given row in Table 1.
a_0_nM = 10
tau_s = np.array([4.16, 3.53])
k = np.array([1.012, .997]) # facilitation > 1 depression < 1
dt_s = 0.01
Vm_uM_sec = 4
Km_uM = .2
I_init = np.array([1, 1]) # If it's 1, then no slow change - seems pointless then.
A_init = .5
C_init = 0.001
#################
n_factors = tau_s.size
tmax_s = 20
t_s = np.arange(start = 0,stop = tmax_s,step = dt_s)
spikes = np.random.rand(t_s.size) > 0.95

I = np.zeros((t_s.size,n_factors))
I[0,] = I_init
C = np.zeros(t_s.size)
A = np.zeros(t_s.size)
C[0] = C_init;
A[0] = A_init;
V = np.zeros(n_factors)

for ii in range(1, t_s.size):
    for jj in range(0,n_factors):
        I[ii,jj] = I[ii-1,jj] + dt_s*(1-I[ii-1,jj])/tau_s[jj]  # here is where I get confused - should it be I - 1 or 1- I - something gets flipped in the conversion the euler method I think

    if spikes[ii]:
        for jj in range(0,n_factors):
            V[jj] = k[jj]*I[ii,jj]
            
        A[ii] = a_0_nM*np.prod(V)

    C[ii] = C[ii-1] + dt_s*(spikes[ii]*A[ii] - Vm_uM_sec/(1 + Km_uM/C[ii-1]))

sns.set_style('white')
plt.figure(figsize=(15,7))
plt.style.use('seaborn-pastel')
plt.plot(t_s, C, color = u'#FFBB6C', label='C')
plt.plot(t_s[spikes], t_s[spikes]*0, '+',color = 'red', label='C')
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('C')


