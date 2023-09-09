#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions including statistics, graphing, and experiment data import
@author: wmarchsteinman
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random
import pandas as pd
import csv
from scipy import interpolate


import os

LV_20 = [0.00, 0.38, 1.01, 1.24]
LV_10 = [0.00, 0.39, 1.00, 1.28]

#Lv via Shinomoto et. al. 2009
def local_variation(I):
    n = len(I)
    sum = 0
    for i in range(0, n-1):
        sum += ((I[i]-I[i+1])/(I[i]+I[i+1]))**2
    return 3/(n-1)*sum

#LvR via Shinomoto et. al. 2009, R = 0.005 refractory constant
def LVR(I, R = 0.005):
    n = len(I)
    sum = 0
    for i in range(0, n-1):
        sum += (1 - (4*I[i]*I[i+1])/(I[i]+I[i+1])**2)*(1 + (4*R)/(I[i]+I[i+1]))
    return 3/(n-1)*sum


#produces spike timings from an ISI list
def generate_timings(I):
    sum = 0
    timings = []
    for i in I:
        sum += i
        timings.append(sum)
        if (sum >= 10):
            break
    return timings

#Generates a list of ISI from a spike timing array
def generate_ISI(T):
    ISI = []
    for i in range (0, len(T)-1):
        ISI.append(T[i+1]-T[i])
    return ISI
    
#Here, I is an ordered list of inter-spike intervals.
def print_stats(I):
    print(f'Local Variation: {local_variation(I)}\nLvR: {LVR(I,.005)}\nMean: {np.mean(I)}\n'
          +f'Variance: {np.var(I)}\nStandard Deviation: {np.std(I)}')


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

#read a timing file into a list
def timing_to_arrayf(filename):
    I = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines:
            I.append(float(l))
    return I

#return local variation from a timing list
def lv_timing(T):
    I = generate_ISI(T)
    return local_variation(I)


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

def trial_average(trials):
    DopAvg = []
    for i in range(4):
        x = np.zeros(len(trials[:,i]))
        for j in range(0,5):
            x += trials[:,i*5 + j]
        DopAvg.append(x/5)
    return DopAvg

    
def graphAverage(data, timings, labels, title, trials = 5, start = 0, stop = 250, error = False, save = False, df = False):
    figure, graph = plt.subplots()
    leng = 0
    if df:
        leng = data.size()
    else:
        leng = len(data)
    for i in range(leng):
        graph.plot(timings[start:stop], data[i][start:stop])
        y_err = 2*np.std(data[i])/np.sqrt(trials)
        if error:
            tar = 4
            time_compressed = []
            t_c = []
            errX_TOT = []
            errY_TOT = []
            for k in range(len(data)):
                time_comp_s = []
                for j in range(len(data[i])):
                    if j%tar == 0:
                        if len(t_c) < len(data[k])/tar:
                            t_c.append(timings[j])
                        time_comp_s.append(data[k][j])
                errX_TOT.append(np.argmax(data[k])/5)
                errY_TOT.append(np.max(data[k]))
                time_compressed.append(time_comp_s)
            graph.errorbar(errX_TOT[i], errY_TOT[i], yerr = y_err, color = 'black', ecolor = 'black', elinewidth = 1, fmt = 'o', capsize = 3, capthick = 1)
        #graph.errorbar(Dop_20_timing[0:250], DopAvg[i][0:250], errorevery = 50, yerr = y_err, elinewidth = 2, capsize = 3, capthick = 1,  )
    graph.set_title(f'{title}')
    graph.set_xlabel('time (s)')
    graph.set_ylabel('Dopamine Concentration (nM)')
    plt.legend(labels)
    plt.show()
    if save:
        plt.savefig('./Graphs_exp/'+ title + '.png', bbox_inches='tight')
        
def graphDataLV(data, LV, title = "Data vs. LV"):
    figure, graph = plt.subplots()
    counter = 0
    for i in data:
        lv = [LV[counter]]*len(i)
        graph.scatter(i, lv, marker = 'x', s = 11)
        counter+=1
        graph.set_yticks(LV)
        graph.set_ylabel('LV')
        graph.set_xlabel('Dopamine Concentration (nM)')
    plt.title(title)
    plt.legend(LV)
    plt.show()
    
        
def graphLVData(data, LV, title = "LV vs. Data", error = False, max = False, tri = 5):
    figure, graph = plt.subplots()
    counter = 0
    for i in data:
        lv = LV[counter]
        d = 0
        y_err = 2*np.std(i)/np.sqrt(tri) #assumes average
        if max == True:
            d = np.max(i)
        else:
            d = np.mean(i)
        if error:
            graph.errorbar(lv, d, yerr = y_err, ecolor = 'black', elinewidth = 1, fmt = 'o', capsize = 3, capthick = 1)
        else:
            lv = [LV[counter]]*len(i)
            graph.scatter(lv, i,marker = 'o', s = 11)
        counter+=1
        graph.set_xticks(LV)
        graph.set_xlabel('LV')
        graph.set_ylabel('Dopamine Concentration (nM)')
    plt.title(title)
    plt.legend(LV)
    plt.show() 

def plotCDF(Data, LV, title = ""):
    figure, graph = plt.subplots()
    for i in range(len(Data)):
        x = np.sort(Data[i])
        y = 1. * np.arange(len(Data[i])) / (len(Data[i]) - 1)
        plt.plot(x, y)
    plt.title(title)
    plt.set_xlabel('Dopamine Concentration (nM)')
    plt.legend(LV)
    plt.show()
    
def plotPDF(Data, LV, title = "", bin = 50):
    figure, graph = plt.subplots()
    for i in range(len(Data)):
        figure, graph = plt.subplots()
        graph.hist(Data[i], density = True, bins = bin)
        plt.title(title)
        graph.set_xlabel('Dopamine Concentration (nM)')
        plt.show()
        

#interpolate 
def interpolate_trace(timings, data, deriv = 0):
    x =np.linspace(timings[0], timings[-1], len(timings)*4)
    f = interpolate.splrep(timings, data, s=0) #interp1d(timings, data, kind='nearest')
    y = interpolate.splev(x, f, der=deriv)
    return x, y

def forward_difference(trace):
    arr = []
    time = [0.1]
    dt = 0.2
    for i in range(1,len(trace)):
        arr.append((trace[i]-trace[i-1])/dt)
        time.append(time[i-1]+dt)
    return time[0:-1], arr

def second_central_difference(trace):
    arr = []
    time = [0.1]
    dt = 0.2
    for i in range(1,len(trace)-1):
        arr.append((trace[i+1]-2*trace[i] + trace[i-1])/dt**2)
        time.append(time[i-1]+dt)
    return time[0:-1], arr

def graph_differences(trace, timing, samps = [], gsamps = False, title = "",smoothed = False, error = False, accel = True, vel = True, main = True, wind = 20, tstart = 0, tend = 500, legend = LV_20):
    for i in range(len(trace)):
        maxPos = np.argmax(trace[i])
        max0 = np.max(trace[i])
        time, forward_diffArr = interpolate_trace(timing, trace[i], deriv = 1)
        time2, scd_arr = interpolate_trace(timing, trace[i], deriv = 2)
        if smoothed:
            df= dv = da = None
            if main:
                df = pd.DataFrame(np.transpose(trace), index = timing).ewm(span = wind).mean()
                df.plot()
            if vel:
                dv = pd.DataFrame(np.transpose(forward_diffArr), index = time).ewm(span = wind).mean()
                dv.plot()
            if accel:
                da = pd.DataFrame(np.transpose(scd_arr), index = time2).ewm(span = wind).mean()
                da.plot()
            if gsamps:
                plt.eventplot(10+np.array(samps)[i], orientation="horizontal", linewidth=0.5, linelengths = 100, lineoffsets = 0, color = "gray")
            #plt.show()
            plt.title(f"LV ={legend[i]}")
            plt.xlabel('time (s)')
            plt.ylabel('Dopamine (nM/s)')
            #plt.clf()
        else:
            figure, graph = plt.subplots()
            if main:
                graph.plot(timing[0:maxPos], trace[i, 0:maxPos])
            if vel:
                graph.plot(time, forward_diffArr)
            if accel:
                graph.plot(time2, scd_arr)
            plt.title(title)
            graph.set_xlabel('time (s)')
            plt.show()




