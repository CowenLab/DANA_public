from matplotlib import pyplot as plt
import numpy as np
import random

# Michaelis–Menten alone - nothing else
tspan_s = 10
dt = 0.001  # time step for forward euler method
loop  = np.ceil(tspan_s/dt)   # no. of iterations of euler
loop = loop.astype(int)

t = np.arange(1,loop)*dt
stim_train = t*0
C = t*0

r = np.unique(np.ceil(np.random.rand(100, 1)*(t.size-100))) # random stim pattern.
stim_train[r.astype(int)] = 1 # maybe this needs to be some other value.

# I have no idea what reasonable values of these are...
Cp = 20 # concentration of DA released with each pulse.
Vm = 5 # maximal velocity
Km = 2 # affinity constant
for i in np.arange(1,loop-2):
    # dCdt = r*Cp - Vm/(1+Km/C);
    C[i+1] = C[i] + Cp*stim_train[i] - Vm/(1+Km/C[i])

plt.figure()
plt.plot(t,C)
plt.xlabel('sec')
plt.ylabel('[DA]')
plt.show()