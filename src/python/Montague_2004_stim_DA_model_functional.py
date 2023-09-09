"""
Montague model. 
Montague PR, McClure SM, Baldwin PR, Phillips PEM, Budygin E a, Stuber GD, et al. Dynamic gain control of dopamine delivery in freely moving animals. The Journal of Neuroscience : The Official Journal of the Society for Neuroscience. 2004;24:1754–1759.

Also added a simple synaptic release model from Miller, Paul. An Introductory Course in Computational Neuroscience

@author: Stephen Cowen
@year: 2022
"""
import numpy as np
import matplotlib.pyplot as plt
import os
#############################################################
# NOTE: Users other than Stephen should only have to modify the data_dir below to point to their particular GitHub path to the DANA folder.
# This is the only way you will be able to load the stim and DA files.
data_dir = r'C:\Users\Stephen Cowen\Documents\GitHub\DANA\Data\Acute\20220218'
stim_fname = 'stim_times_10Hz_LV1pt0.txt' # Assumes a simple text file of just timestamps.
da_fname = 'CV_10Hz_1pt0_trial_1_DOS.csv' # Need to export from excel as DOS .csv (not other types)
######################################################

def Montague_model(a_0_nM = 10,tau_s = np.array([4.16, 3.53]), k = np.array([1.012, .997]), 
                   dt_s = 0.001, tmax_s = 20, Vm_uM_sec = 4, Km_uM = .2, spikes_sec = None,
                   I_init = np.array([1.5, 1.4]), A_init = .5, C_init = 0.001, PLOT_IT = True ):
    import numpy as np
    
    # Use the values from a given row in Table 1.
    #a_0_nM = 10
    #tau_s = np.array([4.16, 3.53]) # FOr some reason, changing this paramter does not seem to affect much - it should REALLY affect things!!!!
    #k = np.array([1.012, .997]) # facilitation > 1 depression < 1
    #dt_s = 0.01
    #tmax_s = 20
    #Vm_uM_sec = 4
    #Km_uM = .2
    #I_init = np.array([1.5, 1.4]) # If it's 1, then no slow change - seems pointless then.
    #A_init = .5 # 
    #C_init = 0.001 # 
    #################
    n_factors = tau_s.size
    t_s = np.arange(start = 0,stop = tmax_s,step = dt_s)
    if (spikes_sec is None):
        # Generate random spikes/stims
        n_spikes = 10*tmax_s # 10 Hz
        spikes = np.random.rand(t_s.size) > 1-n_spikes/t_s.size
        spikes_sec = t_s[spikes]
    else:
        # Find the spike times closest to times the user passed in (if they did)
        ix = []
        for ts in spikes_sec:
            ix.append(np.where(t_s >=ts)[0][0]) # Finds the closest match in the t_s array and returns the index in t_s
        spikes = np.zeros(t_s.size) > 0 # Create a boolean array
        spikes[ix] = True
        
    # Initialize variables.
    I = np.zeros((t_s.size,n_factors))
    I[0,] = I_init
    C = np.zeros(t_s.size)
    A = np.zeros(t_s.size)
    C[0] = C_init
    A[0] = A_init
    V = np.zeros(n_factors)
    
    # This is the main loop that simulates dopamine release/reuptake
    for ii in range(1, t_s.size):
        for jj in range(0,n_factors):
            I[ii,jj] = I[ii-1,jj] + dt_s*((1-I[ii-1,jj])/tau_s[jj])  # here is where I get confused - should it be I - 1 or 1- I - something gets flipped in the conversion the euler method I think
    
        if spikes[ii]:
            for jj in range(0,n_factors):
                V[jj] = k[jj]*I[ii,jj]
                
            A[ii] = a_0_nM*np.prod(V)
    
        C[ii] = C[ii-1] + dt_s*(spikes[ii]*A[ii] - Vm_uM_sec/(1 + Km_uM/C[ii-1]))
    
    if PLOT_IT:
        plot_DA_and_stim(t_s,C,spikes,spikes_sec)
        
    return C, t_s, spikes

def Synaptic_release_model(Vth = -0.035, tau_rise = 0.002, tau_decay = 0.025, 
                   dt_s = 0.001, tmax_s = 25, Vrange = 0.005, Gmax = 1e-9, spikes_sec = None, PLOT_IT = True ):
    import numpy as np
    #################
    # From chapter 5 in computational book. Miller, Paul. An Introductory Course in Computational Neuroscience
    # This model really has only 2 free parameters: tau rise and tau decay.
    # Simpler than the Montague model. Does it do as well?
    #
    # This is a model of synaptic conductance G as a function of depolarization of the terminal.
    # Because presynaptic Ca dynamics can be quite slow (which could mirror slow reuptake) and because
    # synaptic conductance should correlate roughly to DA release, it seemed like a reasonable simplification. The output does not have to be interpreted as G, but as DA perhaps, just modeld as a double exponential.
    # The quesiton is, does this simpler model do a decent enough job of modeling the data?
    # (relative to Montague)
    # 
    V_stim = -0.01 # I don't think this will matter much as long as reasonable value above threshold.
    t_s = np.arange(start = 0,stop = tmax_s,step = dt_s)
    if (spikes_sec is None):
        # Generate random spikes/stims
        n_spikes = 10*tmax_s # 10 Hz
        spikes = np.random.rand(t_s.size) > 1-n_spikes/t_s.size
        spikes_sec = t_s[spikes]
    else:
        # Find the spike times closest to times the user passed in (if they did)
        ix = []
        for ts in spikes_sec:
            ix.append(np.where(t_s >=ts)[0][0]) # Finds the closest match in the t_s array and returns the index in t_s
        spikes = np.zeros(t_s.size) > 0 # Create a boolean array
        spikes[ix] = True
        
    # Initialize variables.
    G = np.zeros(t_s.size)
    
    # This is the main loop that simulates dopamine release/reuptake
    for i in range(1, t_s.size):
    
        if spikes[i]:
            V = V_stim
        else:
            V = -0.06 # subthreshold. I am guessing that as long as it's subthresh, it won't matter.
            
        Gsyn_inf = Gmax/(1+np.exp(-(V-Vth)/Vrange)) # signmoid to prevent infinite I believe
        tau_Gsyn = tau_rise + (tau_decay-tau_rise)/(1 + np.exp((V-Vth)/Vrange)) # left or right sigmoid again to squash values into a reasonable range.
        G[i] = G[i-1] + dt_s*(Gsyn_inf-G[i-1])/tau_Gsyn # Euler

    if PLOT_IT:
        plot_DA_and_stim(t_s,G,spikes,spikes_sec)
        
    return G, t_s, spikes


def plot_DA_and_stim(t_s,C,spikes,spikes_sec = None, t_s_DA = None):
    import matplotlib.pyplot as plt
    #import seaborn as sns
    #plt.style.use('seaborn-pastel')
    #plt.figure(figsize=(15,7))
    fig, ax1 = plt.subplots()
    #sns.set_style('white')
    #sns.set(font_scale=2)
    ax1.plot(t_s, C, color = 'black', label='C')
    ax1.plot(t_s[spikes], t_s[spikes]*0, '+',color = 'red', markersize=20)
    if (spikes_sec is not None):
        ax1.plot(spikes_sec, spikes_sec*0, '^',color = 'blue', markersize=5)
    ax1.grid()
    ax1.set_xlabel('Time (sec)')
    ax1.set_ylabel('C')
    if (t_s_DA is not None):
        ax2 = ax1.twinx()
        ax2.plot(t_s_DA[:,0], t_s_DA[:,1])
        plt.xlim(5 , 25)
        

    
# Run Montague on randomly generated stim trains
C, t_s, spikes = Montague_model(a_0_nM = 11,tau_s = np.array([4.16, 3.53]), k = np.array([1.012, .997]),dt_s = 0.01)

# Run Montague on manually generated stim trains
#C, t_s, spikes = Montague_model(a_0_nM = 11,tau_s = np.array([4.16, 3.53]), k = np.array([1.012, .997]),dt_s = 0.001, spikes_sec = np.array([0.2, 1.2, 1.3, 1.34, 2, 3, 4, 4.1, 5]))

# Run Montague on stim trains loaded from a file

stim_file_path = os.path.join(data_dir, stim_fname)
DA_file_path = os.path.join(data_dir, da_fname)

# Load the stim times from a text file.
stim_times_s = np.loadtxt(stim_file_path)
#stim_times_s = stim_times_s[1] # I should not have to do this, but can't seem to get this to work WITHIN the load file.
tmx = stim_times_s[stim_times_s.size-1] 
C, t_s, spikes = Montague_model(a_0_nM = 11,tau_s = np.array([4.16, 3.53]), k = np.array([1.012, .997]),dt_s = 0.01, tmax_s = 22, spikes_sec = stim_times_s)
C, t_s, spikes = Montague_model(a_0_nM = 11,tau_s = np.array([4.16, 3.53]), k = np.array([3.012, .997]),dt_s = 0.01, tmax_s = 22, spikes_sec = stim_times_s)
C, t_s, spikes = Montague_model(a_0_nM = 10,tau_s = np.array([4.16, 3.53]), k = np.array([3.012, .997]),dt_s = 0.01, tmax_s = 22, spikes_sec = stim_times_s, Km_uM = .2)


# Load the DA measurements (and stim times) from a csv file. Assumes one file per trial and the LV and trial ID are indicated.
t_s_DA = np.loadtxt(DA_file_path, delimiter=",")
stim_times_s = np.loadtxt(stim_file_path)
sC, st_s, sspikes = Synaptic_release_model(spikes_sec = stim_times_s, tau_rise = .150, tau_decay = 1.2, tmax_s = 22,dt_s = 0.01)

plot_DA_and_stim(t_s = t_s, C = C, spikes = spikes, t_s_DA = t_s_DA)
plot_DA_and_stim(t_s = st_s, C = sC, spikes = sspikes, t_s_DA = t_s_DA)

# Compute the correlation and R2 for the model and actual output
# They need to be in the same timebase and we should restrict to a reasonable period surrounding the first and last stim pulse.

pred_C_interp = np.interp(t_s_DA[:,0], t_s, C)
# correlate
time_range_s = [10, 22] # reasonable range from which to calculate the correlation...
ix = np.where((t_s_DA[:,0] >= time_range_s[0]) & (t_s_DA[:,0] <= time_range_s[1]))
corr = np.corrcoef(t_s_DA[ix,1],pred_C_interp[ix])
corr = corr[0,1]
coef = np.polyfit(t_s_DA[ix,1][0],pred_C_interp[ix],1)
poly1d_fn = np.poly1d(coef) 

# plot the relationships...
plt.figure(figsize=(15,7))
plt.subplot(1,2,1)
x =t_s_DA[ix,0][0]
y =t_s_DA[ix,1][0]
yz = (y-np.mean(y))/np.std(y)
pred_C_interp_z = (pred_C_interp[ix]-np.mean(pred_C_interp[ix]))/np.std(pred_C_interp[ix])
pred_C_interp_z = pred_C_interp_z
plt.plot(x,yz, 'bo-')
plt.plot(t_s_DA[ix,0][0],pred_C_interp_z,'g+-')
plt.title('r = %1.2f r^2 = %1.2f' % (corr,corr*corr))
plt.ylabel('z')
plt.xlabel('time (s)')
plt.subplot(1,2,2)
plt.plot(y,pred_C_interp[ix], 'yo', y, poly1d_fn(y), '-k')
plt.xlabel('C')
plt.ylabel('predicted C')
plt.title('least squares fit')

