#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 17:05:19 2022

@author: wmarchsteinman
"""
import numpy as np
import matplotlib.pyplot as plt

#import from csv (sheet from excel file)
def experiment_import(file_name, title = '', date = '', header_rows = 3, timing_column = 0, trials = 5, rounds = 4):
    file_array = np.genfromtxt(file_name, delimiter = ',')[header_rows:]
    file_timing = file_array[:,timing_column]
    dop_current = []
    dop_levels = []
    dop_norm20 = []
    dop_norm10 = []
    
    for i in range(rounds):
        offset = i * trials * 4 + 1*(i+1)
        for j in range(trials):
            dop_current.append(file_array[:,offset+j])
            dop_levels.append(file_array[:,offset+trials + j])
            dop_norm20.append(file_array[:,offset + 2*trials + j])
            dop_norm10.append(file_array[:,offset + 3*trials + j])
            
    return np.array(file_timing), np.array(dop_current).T, np.array(dop_levels).T, np.array(dop_norm10).T, np.array(dop_norm20).T


#returns a timing file as an array.
def timing_to_arrayf(filename):
    I = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines:
            I.append(float(l))
    return I

#create an event graph of timing sequence.
def graph_timings(T):
    figure, graph = plt.subplots()
    graph.eventplot(T, orientation="horizontal", linewidth=0.75, linelengths = 0.5, lineoffsets = 0)
    plt.show()

#export timing to file
def export_timing(T, file_name):
    output_file = open('./timing_output/'+file_name+'.txt', "w")
    np.savetxt(output_file, T, fmt = "%.16f")
    output_file.close()
    return True

#export ISI array to file
def export_ISI(I, file_name):
    output_file = open('./ISI_output/'+file_name+'.txt', "w")
    np.savetxt(output_file, I, fmt = "%.16f")
    output_file.close()
    return True