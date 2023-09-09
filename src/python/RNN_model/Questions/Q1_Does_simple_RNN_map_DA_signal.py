# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 19:49:34 2020
This is a 'package'
Runs an autoencoder on the data.

NOTE: What is the outcome measure here? How do we know the autoencoding works better than some other form of autoencoding?

@author: Stephen Cowen 2020.
"""
#%%
import os
#import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import scipy as sp
import glob
import LD_Lib   # Assumes this was added to the sys.path (in LD_Session_Iterator.py)
import LD_ML    # Assumes this was added to the sys.path

PLOT_IT = True


def main(data_dir):
    return(Analysis(data_dir))

def Analysis(data_dir):

    tmp, day = os.path.split(data_dir)
    tmp, animal = os.path.split(tmp)

    os.chdir(data_dir) # go to the data directory

    binsize_ms = 100
    epoch = 1 # zero based indices in python - need to remember this
    OUT = {'data_dir': data_dir, 'binsize_ms': binsize_ms }
    
    ##############################
    # Load the epoch times and position data.
    ##############################
    epochtimes = pd.read_csv("Epochs_uS.csv")   
    trialtimes = pd.read_csv("Trial_Times_uS.csv")   
    posdata = pd.read_csv("Position.csv")
    riptimes = pd.read_csv("Ripples_uS.csv",header=None)
    # Restrict the postion and spike data to the timeframe in question.
    posdata = posdata.loc[(posdata['Time_uS']>epochtimes['Start'][epoch]) & (posdata['Time_uS']<epochtimes['End'][epoch])]
    #spiketime_edges_uS = np.arange(epochtimes['Start'][0], epochtimes['End'][5], binsize_ms*1e3)
    spiketime_edges_uS = np.arange(epochtimes['Start'][epoch], epochtimes['End'][epoch], binsize_ms*1e3)
    
    ##############################
    # Load all of the spike files.
    ##############################
    spk,spk_files = LD_Lib.Load_Spike_Files("NRN*.csv", path = '')
    # Bin the spikes
    NT = LD_Lib.Bin_Spike_Times(spk,spiketime_edges_uS)
    # Determine the meaning of each row in the spike by time matrix (trial type, running speed...)
    # This dataframe is INCREDIBLY important as it has the trial type, trial number, shock or no shock, etc...
    x_binctrs_uS = spiketime_edges_uS[0:-1]+ binsize_ms*1000/2
    df = LD_Lib.Table_From_Times(x_binctrs_uS, epochtimes, posdata, trialtimes, riptimes)
    plt.figure()
    plt.scatter(df['time_uS']/60e6, df['Linear'], c=df["BehaviorState"].cat.codes, s=10, cmap="tab10")
    
    IX = (df["Speed_cmSec"]>= 5).values
    NTrun = NT[IX,:]
    Yrun = df.loc[IX,'LinearByDir'].values
    Yrun = Yrun.reshape(len(Yrun),1)

    ms1,mod1 = LD_ML.Linear_Regression(NTrun,Yrun,PLOT_IT=True)
    ms2,mod2 = LD_ML.Decision_Tree_Regression(NTrun,Yrun,PLOT_IT=True)
    ms3,mod3 = LD_ML.SVM(NTrun,Yrun,kernel = 'linear', PLOT_IT=True)
    ms4,mod4 = LD_ML.SVM(NTrun,Yrun,kernel = 'rbf',PLOT_IT=True)
    ms5,mod5 = LD_ML.MLP(NTrun,Yrun,PLOT_IT=True)
    #ms6,mod6 = LD_ML.Naive_Bayes_Gaussian(NTrun,Yrun,PLOT_IT=True) # Work to do here. need to define the priors.

    from sklearn.model_selection import cross_val_score
    scores = cross_val_score(mod1, NTrun, Yrun,
                         scoring="neg_mean_squared_error", cv=10)
    mod1_scores = np.sqrt(-scores)
    scores = cross_val_score(mod2, NTrun, Yrun,
                         scoring="neg_mean_squared_error", cv=10)
    mod2_scores = np.sqrt(-scores)
    scores = cross_val_score(mod3, NTrun, Yrun,
                         scoring="neg_mean_squared_error", cv=10)
    mod3_scores = np.sqrt(-scores)
    scores = cross_val_score(mod4, NTrun, Yrun,
                         scoring="neg_mean_squared_error", cv=10)
    mod4_scores = np.sqrt(-scores)
    scores = cross_val_score(mod5, NTrun, Yrun,
                         scoring="neg_mean_squared_error", cv=10)
    mod5_scores = np.sqrt(-scores)
    # Plot the data.
    #pt = LD_Lib.Plot_Spike_Time_Matrix(x_binctrs_uS,NT)
    #pt.show()
    #pt.close()
    ###################################################################################################
    # Plot the spike x time matrix for active behavioral episodes and when the rat is actively running.
    # This shows the usefulness of the df: Choose only the data where the rat is running fast and meets all of the other conditions.
    ###################################################################################################
    Indices = {
        "Dir1_No_Stim" : ((df["RunningDirection"]==1) & (df["BadTrial"]==False) & (df["Stimmed"]==False) & (df["Epoch"]=='Maze1') & (df["Speed_cmSec"]>= 5)).values,
        "Dir1_Stim" : ((df["RunningDirection"]==1) & (df["BadTrial"]==False) & (df["Stimmed"]==True) & (df["Epoch"]=='Maze1') & (df["Speed_cmSec"]>= 5)).values,
        "Dir2_No_Stim" : ((df["RunningDirection"]==2) & (df["BadTrial"]==False) & (df["Stimmed"]==False) & (df["Epoch"]=='Maze1') & (df["Speed_cmSec"]>= 5)).values,
        "Dir2_Stim" : ((df["RunningDirection"]==2) & (df["BadTrial"]==False) & (df["Stimmed"]==True) & (df["Epoch"]=='Maze1') & (df["Speed_cmSec"]>= 5)).values
    }

    if PLOT_IT:
        fig = plt.figure()

        for cnt,key in enumerate(Indices):
            plt.subplot(2, 2, cnt+1)
            pt = LD_Lib.Plot_Matrix(NT[Indices[key],:])
            pt.title(key)

        fig.suptitle('Neuron by time matrices when rat is running ' + animal + ' ' + day , fontsize=14)
        pt.show(block=True)
        plt.savefig('C:/Temp/' + animal + day + 'TimeXNeuron.png')

    n_neurons = len(spk)
    OUT.update({'n_neurons': n_neurons})
    n_spikes = np.sum(NT)
    OUT.update({'n_spikes': n_spikes})

    # Do some machine learning for fun...
    PCscores,a,b = LD_ML.PCA(NT, PLOT_IT = True)
    #LD_ML.Kernel_PCA(NT,PLOT_IT = True)
    #c = LD_ML.Kmeans(NT,range(2,20), PLOT_IT = True)
    
    LD_ML.TSNE(NT, df["BehaviorState"].cat.codes, PLOT_IT = True)

    LD_ML.TSNE(PCscores[:,1:10], df["BehaviorState"].cat.codes, PLOT_IT = True)
    #aNT = LD_ML.Autoencode_ST_Matrix_1Layer(NT)
    print('Stacked')
    NTcoded = LD_ML.Autoencode_Stacked(NT, n_units = (np.round(n_neurons/2),8,4), epochs = 20, PLOT_IT = True)
    #aNT = LD_ML.Autoencode_ST_Matrix_Stacked(PCscores[:,0:6],PLOT_IT = True)
    LD_ML.TSNE(NTcoded, df["BehaviorState"].cat.codes, PLOT_IT = True)
    plt.show()

    return(OUT)

def Meta_Analysis(meta_dir):
    # This function is here to analyze all of the output of the Analysis function and produce summary plots
    #
    os.chdir(meta_dir)
    print("Starting Meta Analysis. Loading .pkl files.")
    data_files = glob.glob(meta_dir + "/*.pkl")
    for data_file in data_files:
        print (data_file)
        df = pd.read_pickle(data_file)
        print(df)


if __name__ == "__main__":
    main(r'C:\Temp\OK\8417\Day10') # main(os.getcwd())
    #main(int(sys.argv[1]))

# %%
