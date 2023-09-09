#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 16:52:45 2022

@author: wmarchsteinman
"""
import numpy as np
import matplotlib.pyplot as plt

from pandas import read_csv
from matplotlib import pyplot
from pandas.plotting import lag_plot
from pandas import DataFrame, concat
from pandas.plotting import autocorrelation_plot
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.ar_model import AutoReg
from statsmodels.tsa.arima.model import ARIMA
from sklearn.metrics import mean_squared_error
from math import sqrt
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


import os
def timing_to_arrayf(filename):
    I = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines:
            I.append(float(l))
    return I

file_list = os.listdir('./timing_data/timings_used/')

samples = []
counter = 0
samples_order = [.38,1.01, 1.24, .39, 0, 1.00, 1.28]
demo = np.sort(np.random.rand(11)*10)
print_stats(demo)
graph_timings(demo)

demo_0 = np.linspace(0, 10, 11)
print_stats(demo_0)
graph_timings(demo_0)

for name in file_list:
    timings = timing_to_arrayf('./timing_data/timings_used/'+name)
    samples.append(timings)
    intervals = generate_ISI(timings)
    #print_stats(reorder(intervals))
    print(f'index: {counter}, Name: {name}')
    #print_stats(intervals[:int(len(intervals)/2)])
    #print_stats(intervals[int(len(intervals)/2):])
    #graph_timings(timings)
    
    counter+=1


spikes10 = [samples[4], samples[3], samples[5], samples[6]]
spikes20 = [np.linspace(0, 10, 200), samples[0], samples[1], samples[2]]


T_R = 2
V_m = 4
K_m = .2
a_0 = 10 #irregular interburst interval from Montague 2004
c = [.001]
coeffs = [.45, 1.01, .35, .984, 600, .99] #standard values from paper -- T1 k1 T2 k2 T3 k3 -- can add more if necessary in pairs

results = 0
solutions_LV20 = []
end = 20
#spikes20[1] = np.zeros(len(spikes20[1]))
#spikes20[1][0] = 1
for i in range(0,len(spikes20)):
    spike_timings = spikes20[i]
    res = solve_ivp(Montague_model, (0,end), c, method = 'RK45', 
                  args = (V_m, K_m, a_0, coeffs, spike_timings, T_R), dense_output = True)
    solutions_LV20.append(res)

solutions_LV10 = []
c = [0.001]
for i in range(len(spikes10)):
    spike_timings = spikes10[i]
    res = solve_ivp(Montague_model, (0,end), c, method = 'RK45', 
                  args = (V_m, K_m, a_0, coeffs, spike_timings, T_R), dense_output = True)
    solutions_LV10.append(res)
    
t = np.linspace(0, 20, 300)
figure, graph = plt.subplots()
#graph.plot(t, sol.sol(t)[0])
for i in solutions_LV20:
    graph.plot(t, i.sol(t)[0]*80)

#testing versus 4/1/22 20Hz normalized
plt.gca().set_prop_cycle(None)
Avg_Norm_2020 = trial_average(norm2020_31822)
time_domain = np.linspace(0, 20, 101)

for i in Avg_Norm_2020:
    graph.plot(time_domain, i[50:151],linestyle = 'dashed')

plt.ylabel('Dopamine Concentration (uM)')
plt.xlabel('time (s)')
plt.title('Modeled Dopamine release (Montague) for 20Hz spike train timings')
plt.legend(LV_20)
plt.show()

figure, graph = plt.subplots()
for i in solutions_LV10:
    graph.plot(t, i.sol(t)[0]*550)

plt.gca().set_prop_cycle(None)
Avg_Norm_1010 = trial_average(norm1010_4822)
time_domain = np.linspace(0, 20, 101)
for i in Avg_Norm_1010:
    graph.plot(time_domain, i[50:151],linestyle = 'dashed')
    
plt.ylabel('Dopamine Concentration (uM)')
plt.xlabel('time (s)')
plt.title('Modeled Dopamine release (Montague) for 10Hz spike train timings')
plt.legend(LV_10)
plt.show()



