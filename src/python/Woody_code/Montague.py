#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Montague Model
Montague et. al., Dynamic Gain Control of Dopamine Delivery in Freely Moving Animals, 2004

@author: wmarchsteinman
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#mean spike dopamine release -- grant
def mean_spike_dop_release(t, k):
    p_0 = 0.01
    N_k = .001
    gamma = 1000
    n_s = 1
    for i in range(0, k):
        N_k += (1 - p_0)*N_k * np.exp(-gamma*(t[i+1] - t[i])) + p_0*n_s*(1 - np.exp(-gamma*(t[i+1] - t[i]))) #- np.exp(-gamma*(t[i+1] - t[i])/2)*N_k
    return N_k

#total dopamine prediction via grant application
def total_d_k(t, k):
    tot = 0
    p_0 = 0.01
    N_k = .001
    gamma = 1000
    n_s = 1
    for i in range(0, k):
        tot += N_k - tot*np.exp(0.1*(t[i+1] - t[i]))
        N_k += (1 - p_0)*N_k * np.exp(-gamma*(t[i+1] - t[i])) + p_0*n_s*(1 - np.exp(-gamma*(t[i+1] - t[i]))) #- np.exp(-gamma*(t[i+1] - t[i])/2)*N_k
    return N_k, tot

def Wightman_model(x, r, C_p, V_m, K_m):
    return r*C_p - V_m*(1/(1 + K_m/x[0])) - x[1]

#here, we have the following arguments:
#x is a vector containing current dopamine level
#spike_timings is a list of the pulse initiation times
#coeffs = [T_1, k_1, T_2, k_2, T_3, k_3]
#T_i is time constant associated with ith exponential
#k_i is kick associated with ith exponential, k_i > 1 -> facilitation, k_i < 1 -> depression
#V_m is the maximal velocity for dopamine reuptake via Michaelis-Menten kinetics
#empirically-determined value for V_m is 4.0um/sec
#K_m is the associated affinity constant
#empirically-determined value for K_m is 0.2um
#T_R is the time constant of the electrochemical probe used (low-pass filtering via Bath et al., 2000).
#output is concentration in uM
def Montague_model(t, x, spike_timings, coeffs, V_m = 4, K_m = .2, a_0 = 10, T_R = 2):
    k = (p(t, spike_timings, T_R)*A(t, a_0, coeffs) - V_m*(1/(1 + K_m/x[0])))
    #if (t > 10):
        #print(f't: {t}, k: {k}, x: {x[0]}, p: {p(t, spike_timings, T_R)}, A: {A(t, a_0, coeffs)}, second_term: {V_m*(1/(1 + K_m/x[0]))}')
    return k

def p(t, spike_timings, T_R):
    sum = 0
    for i in spike_timings:
        #pay attention only to spikes in the past.
        if (t - i) >= 0:
            sum += eps(t, T_R, i)
    return sum
    #models exact impulse times

#A function via Montague et. al 2004
#Dynamic gain control of release, limited here to three exponential "hidden" dynamic processes
def A(t, a_0, coeffs):
    I_list = []
    for i in range(0,int(len(coeffs)/2)):
        T = coeffs[i*2]
        k = coeffs[i*2 + 1]
        I_list.append(k*(1 - np.exp((-t/T))))
    
    return a_0*np.prod(I_list)
    
def eps(t, T_R, t_i):
    return 1/T_R*np.exp(-(t - t_i)/T_R)


def solveMont(params, spike_timings, t_int = (0, 20), c = 0.001):
    a = solve_ivp(Montague_model, t_int, [c], method = 'RK45', 
                  args = (spike_timings, params[4:], params[1], params[2], params[3], params[0]), 
                  dense_output = True)
    return a

def solveMontDefault(params, spike_timings, t_int = (0, 20), c = 0.001):
    a = solve_ivp(Montague_model, t_int, [c], method = 'RK45', 
                  args = (spike_timings, params), 
                  dense_output = True)
    return a

def graph_sample_mean_spike_dop(sample):
    dop_vals = []
    for i in range(0, len(sample)):
        #dop_vals.append(mean_spike_dop_release(sample, i))
        dop_vals.append(total_d_k(sample, i)[1])
    figure, graph = plt.subplots()
    graph.plot(sample, dop_vals)
    plt.show()
