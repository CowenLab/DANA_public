#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some examples of library usage, including timing file import and Montague model usage

@author: wmarchsteinman
"""
import numpy as np
import matplotlib.pyplot as plt
import os

from scipy.integrate import solve_ivp

from ExperimentFileImport import *
from Utilities import *
from Montague import *
from Visualizations import *

#Data imports

#file imports, of the form: timing, current, dopamine level, first normalization, second normalization (re experiment excel files)
t10_4822, current10_4822, dop10_4822, norm1010_4822, norm1020_4822 = experiment_import('./Experiment_files_4-18/20220408_10Hz.csv')
t20_4822, current20_4822, dop20_4822, norm2010_4822, norm2020_4822 = experiment_import('./Experiment_files_4-18/20220408_20Hz.csv')
t10_4122, current10_4122, dop10_4122, norm1010_4122, norm1020_4122 = experiment_import('./Experiment_files_4-18/20220401_10Hz.csv')
t20_4122, current20_4122, dop20_4122, norm2010_4122, norm2020_4122 = experiment_import('./Experiment_files_4-18/20220401_20Hz.csv')
t10_31822, current10_31822, dop10_31822, norm1010_31822, norm1020_31822 = experiment_import('./Experiment_files_4-18/20220318_10Hz.csv')
t20_31822, current20_31822, dop20_31822, norm2010_31822, norm2020_31822 = experiment_import('./Experiment_files_4-18/20220318_20Hz.csv')
t10_32822, current10_32822, dop10_32822, norm1010_32822, norm1020_32822 = experiment_import('./Experiment_files_4-18/20220328_10Hz.csv')
t20_32822, current20_32822, dop20_32822, norm2010_32822, norm2020_32822 = experiment_import('./Experiment_files_4-18/20220328_20Hz.csv')




file_list = os.listdir('./timing_data/timings_used/') #directory for timing files

samples = []

#read spike timings from files
for name in file_list:
    timings = timing_to_arrayf('./timing_data/timings_used/'+name)
    samples.append(timings)
    intervals = generate_ISI(timings)

#storing spike timings in order of LV for 10Hz and 20Hz trials
spikes10 = [samples[4], samples[3], samples[5], samples[6]]
spikes20 = [np.linspace(0, 10, 200), samples[0], samples[1], samples[2]]


T_R = 2
V_m = 4
K_m = .2
a_0 = 10 #irregular interburst interval from Montague 2004
c = [.001] #initial concentration
coeffs = [4.5, 1.01, 3.5, .984, 600, .99] #standard values from paper -- T1 k1 T2 k2 T3 k3 -- can add more if necessary in pairs

solutions_LV20 = []
t_start = 0
t_end = 20
#spikes20[1] = np.zeros(len(spikes20[1]))
#spikes20[1][0] = 1
for i in range(0,len(spikes20)):
    spike_timings = spikes20[i]
    res = solve_ivp(Montague_model, (t_start,t_end), c, method = 'RK45', 
                  args = (spike_timings, coeffs, V_m, K_m, a_0, T_R), dense_output = True)
    solutions_LV20.append(res)

solutions_LV10 = []
c = [0.001]
for i in range(len(spikes10)):
    spike_timings = spikes10[i]
    res = solve_ivp(Montague_model, (t_start,t_end), c, method = 'RK45', 
                  args = (spike_timings, coeffs, V_m, K_m, a_0, T_R), dense_output = True)
    solutions_LV10.append(res)
    
    
    
#plotting solutions
t = np.linspace(0, 20, 300)
figure, graph = plt.subplots()
for i in solutions_LV20:
    graph.plot(t, i.sol(t)[0]*1000) #scale from uM to nM, note that we use the first list element


plt.gca().set_prop_cycle(None)
Avg_Norm_2020 = trial_average(dop20_31822)
time_domain = np.linspace(0, 20, 101)

for i in Avg_Norm_2020:
    graph.plot(time_domain, i[50:151],linestyle = 'dashed')

plt.ylabel('Dopamine Concentration (nM)')
plt.xlabel('time (s)')
plt.title('Modeled Dopamine release (Montague) for 20Hz spike train timings')
plt.legend(LV_20)
plt.show()

figure, graph = plt.subplots()
for i in solutions_LV10:
    graph.plot(t, i.sol(t)[0]*1000)

plt.gca().set_prop_cycle(None)
Avg_Norm_1010 = trial_average(dop10_31822)
time_domain = np.linspace(0, 20, 101)
for i in Avg_Norm_1010:
    graph.plot(time_domain, i[50:151],linestyle = 'dashed')
    
plt.ylabel('Dopamine Concentration (nM)')
plt.xlabel('time (s)')
plt.title('Modeled Dopamine release (Montague) for 10Hz spike train timings')
plt.legend(LV_10)
plt.show()


#for demonstration only, not valid spike timing lists
#generate random spike timing list, 10Hz avg
demo = np.sort(np.random.rand(100)*10)
print_stats(generate_ISI(demo))#printing stats of an ISI list
graph_timings(demo)

#generate Lv 0 spike timing list, 10Hz avg
demo_0 = np.linspace(0, 10, 100)
print_stats(generate_ISI(demo_0))
graph_timings(demo_0)

#example of using built in functions

#using solveMontDefault()
t_int = (0, 20) #modeling interval
c_0 = .001 #starting dopamine, in uM
#coeffs from above, uses default values for other parameters
result = solveMontDefault(coeffs, demo_0, t_int, c_0)
figure, graph = plt.subplots()
t = np.linspace(0, 20, 300) #evaluate result at at these points -- result is continuous
graph.plot(t, result.sol(t)[0])
plt.ylabel('Dopamine Concentration (uM)')
plt.xlabel('time (s)')
plt.show()


#example of specifying other parameters in solveMont()
T_R = 2
V_m = 4
K_m = .2
a_0 = 10 #irregular interburst interval from Montague 2004
c = [.001] #initial concentration
coeffs = [4.5, 1.01, 3.5, .984, 600, .99]

params = [T_R, V_m, K_m, a_0, 4.5, 1.01, 3.5, .984, 600, .99]

t_int = (0, 20) #modeling interval
c_0 = .001 #starting dopamine, in uM
#coeffs from above, uses default values for other parameters
result = solveMont(params, demo_0, t_int, c_0)
figure, graph = plt.subplots()
t = np.linspace(0, 20, 300) #evaluate result at at these points -- result is continuous
graph.plot(t, result.sol(t)[0])
plt.ylabel('Dopamine Concentration (uM)')
plt.xlabel('time (s)')
plt.show()
