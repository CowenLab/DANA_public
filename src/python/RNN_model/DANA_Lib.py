# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 19:49:34 2020

Collection of functions for the neural ensemble analyis project.

@author: Stephen Cowen 2020.
"""
#%%
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Global parameters for all sessions (e.g., color for plots, ?)
GP = {}

def Load_Spike_Files(identifier_string, path = ''):

    spk = []
    spk_files = glob.glob(path + identifier_string)
    for fn in spk_files:
        tmp = pd.read_csv(fn,header=None).values
        spk.append(tmp)

    return(spk, spk_files)

def Bin_And_Smooth_Spike_Times(all_neurons_us, start_end_us, small_bin_us, big_bin_us, window_size_us):
    #
    from scipy.interpolate import interp1d # https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
    import scipy.signal 
    PLOT_IT = False
    small_spiketime_edges_uS = np.arange(start_end_us[0], start_end_us[1], small_bin_us)
    small_spiketime_ctrs_uS = small_spiketime_edges_uS[0:(len(small_spiketime_edges_uS))-1] + small_bin_us/2

    big_bin_centers_uS = np.arange(start_end_us[0]+big_bin_us/2, start_end_us[1]-big_bin_us/2, big_bin_us)

    if window_size_us is None:
        window_size_us = big_bin_us*2

    #big_bin_n = round(big_bin_us/small_bin_us)

    NT = np.zeros([len(big_bin_centers_uS), len(all_neurons_us)], dtype = 'float16') # Neuron x time matrix. Neuron = column (guess it should be a TN not NT but whatever)

    knl = np.hamming(int(2*window_size_us/small_bin_us)+1) # Decay
    #knl = np.ones(int(2*window_size_us/small_bin_us)+1) # Moving average - half
    
    mid = round(len(knl)/2+.0001)
    knl[:int(mid-1)] = 0
    knl = knl/np.sum(knl)

    #H,ed = np.histogram(all_neurons_us[0], bins=(small_spiketime_edges_uS))
    # Convert to instantaneous firing rates. This makes interpretation a little easier.
    #H = H/(small_bin_us/1e6)
    # Convolve at the small window
    #Hs = np.convolve( H, knl, mode = 'same' )

    # Perform a moving mean on this - but use a half window as we don't want to have spikes providing information about the past as that would impute the spike with the capacity to predict things 1/2 binwidth before the spike even occurred..
    #mid = round(big_bin_n/2+.0001)

    #knl = np.ones((big_bin_n,))
    #knl[:int(mid-1)] = 0
    #knl = knl/sum(knl)
    #Hss = np.convolve(Hs, knl, mode='same')
    

    #fx = interp1d(small_spiketime_ctrs_uS, Hss, bounds_error = False)
    #Hb = fx(big_bin_centers_uS)

    for count,spks in enumerate(all_neurons_us):
        H,ed = np.histogram(spks, bins=(small_spiketime_edges_uS))
        # Convert to instantaneous firing rates. This makes interpretation a little easier.
        H = H/(small_bin_us/1e6)
        Hs = np.convolve( H, knl, mode = 'same' )
        fx = interp1d(small_spiketime_ctrs_uS, Hs, bounds_error = False)
        NT[:,count] = fx(big_bin_centers_uS)


    if PLOT_IT:
        plt.figure()
        plt.plot(small_spiketime_ctrs_uS, H,'.-')
        plt.plot(small_spiketime_ctrs_uS, Hs,'.-')
        plt.plot(big_bin_centers_uS, NT[:,count],'.-')
        plt.figure()
        Plot_Matrix(NT[:1000,:])
        plt.show()


    return (NT,big_bin_centers_uS)

def Bin_Spike_Times(all_neurons,spiketime_edges_uS):
    NT = np.zeros([spiketime_edges_uS.size - 1,len(all_neurons)],dtype = 'float16') # Neuron x time matrix. Neuron = column (guess it should be a TN not NT but whatever)

    for count,spks in enumerate(all_neurons):
         # because path is object not string
         H,edges = np.histogram(spks, bins=(spiketime_edges_uS))
         NT[:,count] = H

    return (NT)

def Plot_Spike_Time_Matrix(t_uS,NT):
    colormap = 'jet'
    plt.cla()
    plt.pcolor(X = t_uS/1e6, C = NT.T, cmap = colormap) 
    #plt.pcolor(X = t_uS/1e6,np.arange(0,NT.shape[0]),NT, cmap = colormap) 
    #plt.imshow(NT, extent = (0,1,0,1), cmap = colormap) # Super fast compared to pcolor but inflexible
    plt.colorbar()
    plt.ylabel('Neuron ID')
    plt.xlabel('Time (s)')    
    return (plt)

def Plot_Matrix(NT):
    colormap = 'jet'
    plt.cla()
    plt.pcolor(NT.T, cmap = colormap) 
    #plt.imshow(NT, extent = (0,1,0,1), cmap = colormap) # Super fast compared to pcolor but inflexible
    plt.colorbar()
    plt.ylabel('Neuron ID')
    plt.xlabel('x')    
    return (plt)

def Plot_NT_Side_By_Side(NT1,NT2):
    colormap = 'jet'
    plt.subplot(1,2,1)
    plt.pcolor(NT1, cmap = colormap) 
    plt.xlabel('Neuron ID')
    plt.subplot(1,2,2)
    plt.pcolor(NT2, cmap = colormap) 
    plt.xlabel('Neuron ID')
    plt.colorbar()



def Table_From_Times(times_uS, epochtimes, posdata, trialtimes, riptimes):
    import pandas as pd
    from scipy.interpolate import interp1d # https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html

    # First let's put all of the data into a dictionary. Maybe this is backwards but seems 
    df = pd.DataFrame()
    df['time_uS'] = times_uS
    fx = interp1d(posdata['Time_uS'],posdata['X'],bounds_error = False) # boundserror allows it to still work if x goes beyond f()
    fy = interp1d(posdata['Time_uS'],posdata['Y'],bounds_error = False)
    fs = interp1d(posdata['Time_uS'],posdata['Speed_cmSec'],bounds_error = False)
    flin = interp1d(posdata['Time_uS'],posdata['Linear'],bounds_error = False)
    df['X'] = fx(df['time_uS'] )
    df['Y'] = fx(df['time_uS'] )
    df['Speed_cmSec'] = fs(df['time_uS'] )
    df['Linear'] = flin(df['time_uS'] )
    df = df.astype({"X": 'float16',"Y": 'float16',"Speed_cmSec": 'float16',"Linear": 'float16'})

    # Go through each row in epoch times and assign a name to that row in the data. df["A"].astype('category') dfc.loc[0, 'A'] = 11
    df['Epoch'] = ''
    for i in range(0,df.shape[1]):
        #print([epochtimes['Epoch'][i],epochtimes['Start'][i], epochtimes['End'][i] ])
        df.loc[(df['time_uS']  > epochtimes['Start'][i]) & (df['time_uS'] < epochtimes['End'][i]), 'Epoch'] = epochtimes['Epoch'][i]

    df = df.astype({"Epoch": 'category'})

    vbls = ['TrlNum', 'RunningDirection','Stimmed', 'BadTrial' ]
    
    #np.dtype(trialtimes.RunningDirection,'int8')
    for iVbl in range(0,len(vbls)):
        for iRow in range(0,trialtimes.shape[0]):
            #print([vbls[iVbl], iRow])
            df.loc[(df['time_uS']  > trialtimes['TrlStartEnd_1'][iRow]) & (df['time_uS'] < trialtimes['TrlStartEnd_2'][iRow]), vbls[iVbl]] = trialtimes[vbls[iVbl]][iRow]


    # Create some meaningful aggregate categories...
    df['BehaviorState'] = ''
    Indices = {
        "Dir1_No_Stim" : ((df["RunningDirection"]==1) & (df["BadTrial"]==False) & (df["Stimmed"]==False) & (df["Speed_cmSec"]>= 5)).values,
        "Dir1_Stim" : ((df["RunningDirection"]==1) & (df["BadTrial"]==False) & (df["Stimmed"]==True) & (df["Speed_cmSec"]>= 5)).values,
        "Dir1_Slow" : ((df["RunningDirection"]==1) & (df["BadTrial"]==False) & (df["Speed_cmSec"] < 5)).values,
        "Dir2_No_Stim" : ((df["RunningDirection"]==2) & (df["BadTrial"]==False) & (df["Stimmed"]==False) & (df["Speed_cmSec"]>= 5)).values,
        "Dir2_Stim" : ((df["RunningDirection"]==2) & (df["BadTrial"]==False) & (df["Stimmed"]==True) & (df["Speed_cmSec"]>= 5)).values,
        "Dir2_Slow" : ((df["RunningDirection"]==2) & (df["BadTrial"]==False) & (df["Speed_cmSec"] < 5)).values,
        "Rew_Zone1" : ((df['Linear']>0) & (df['Linear']<GP['RewardZoneLinPos'][0]) & (df["Speed_cmSec"] < 5)).values,
        "Rew_Zone2" : ((df['Linear']>GP['RewardZoneLinPos'][1]) & (df['Linear']<360) & (df["Speed_cmSec"] < 5)).values
    }

    for key in Indices:
        df.loc[Indices[key],'BehaviorState'] = key

    df = df.astype({"TrlNum": 'int16', "RunningDirection": 'int8', "Stimmed": bool, "BadTrial": bool,"BehaviorState": 'category'}) # some data has infs or NA in it which causes to crash - need to figure out why

    df.loc[:,'BehaviorCodes'] = df["BehaviorState"].cat.codes
    # Create a new variable that has a unique location for each running direction. 
    df.loc[:,'LinearByDir'] = df.loc[:,('Linear')]
    df.loc[df["RunningDirection"]==2,'LinearByDir'] = df.loc[df["RunningDirection"]==2,('Linear')]*-1 + 720

    # plt.plot(df['LinearByDir'],'o')
    # plt.plot(df['Linear'],'r.')
    # plt.show()
    # Bin locations on the track.

    return(df)

def Create_Training_Test_Set_By_Trial(df, prop_in_test_set = 0.5):
    import pandas as pd
    # in order to preserve information from the time-series of activity, we need to pull out contihous periods of time within a trial.
    # I think that the best storage would be a list of lists of indices. Could also do a dictionary. a = ((2,1,2),(2,1))
    n_trials = max(df['Trial'])
    shuffled_indices = np.random.permutation(n_trials)
    test_set_size = int(n_trials * prop_in_test_set)
    test_indices = shuffled_indices[:test_set_size]
    train_indices = shuffled_indices[test_set_size:]
    for iTr in range(0,n_trials):
        print(iTr)

    return(train_indices, test_indices)