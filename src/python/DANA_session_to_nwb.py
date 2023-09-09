#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
See DANA_pynwb_demo.py and run it if it's your first time. 
If this demo does not work, then this code will not work.

The ultimate goal is to upload the data to DANDI
https://www.dandiarchive.org/handbook/10_using_dandi/#uploading-a-dandiset

Our dandi archive/identifier is 
"Brain_stim_and_FSCV_and_ensemble_recording_anesthetized"
I created it on July 28 2022
https://dandiarchive.org/dandiset/000298

To get things to upload to DANDI, you need a key. This is my key (from my DANDI account)
b13214790dba2cbed7cdb5d48d429920d39778a6

Here is another one from staging.
76b5a654c9f1c329fee53bd197a1ce5712e974e1

God damn - between NWB and DANDI - it's a bunch of CS majors with zero conception of how much time their 
impendetrable/buggy software and hours-long video tutorials are stealing from faculty and students. 

For testing, upload and work with datasets here https://gui-staging.dandiarchive.org/


@author: Cowen

"""
# Required libraries. If it crashes here, it means you did not install the library (using conda or pip)
# Be sure to log in as admin to anaconda for installing most packages. Need admin privilidges.
from datetime import datetime
from dateutil.tz import tzlocal
import numpy as np
# The following are specific to nwb
from pynwb import NWBFile, NWBHDF5IO
from pynwb.ecephys import ElectricalSeries, LFP, TimeSeries
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
#
import matplotlib.pyplot as plt # for plotting
from nwbwidgets import nwb2widget 
##########################################################################
# THE FOLLOWING ARE VALUES THAT THE EXPERIMENTER MUST TAILOR TO THE DATASET:
##########################################################################
#
# Decide what data directory we are going to analyze.
session_data_dir = r'Z:\Data\ACUTE_DANA\220408_MFB_NAc'
# For possible fields you can use, see https://pynwb.readthedocs.io/en/stable/pynwb.file.html
nwbfile = NWBFile(
    experiment_description='Determine if different forms of inter-stim-pulse variability alters the magnitude of dopamine release',
    notes='animal did not go under anesthesia for about 30 min. Weak FSCV evoked release. We had about 10 neurons.',
    session_description='FSCV and single unit anesthetized with MFB stim all LV value',
    identifier='220408_MFB_NAc', # For now, can be same as session ID. It is a unique for this nwb file (in theory, there could be mulitple nwb files for a single session (e.g., one for LFP, another for FSCV)
    stimulus_notes ='Stimulation trains delivered at 300uA to the MFB', # 
    session_start_time = datetime(2022, 4, 8, tzinfo=tzlocal()), # Time approximately when the recording actually began
    session_id='220408_MFB_NAc', # (optional) Unique identifier for this session/rat. The session folder seems like a good UID.
    experimenter='Abhi Vishwanath', # (optional)it only allows you to enter ONE experimenter. Ugh. Put one of the experimenters here.
    lab='Cowen/Heien Laboratories', # (optional)
    institution='University of Arizona' # (optional)
)
# Add subject information.
nwbfile.subject = Subject(    
    age='3mo', # for the rat, not the experimenter.
    species='rat', # probably need the propoer scientific name but whatever.
    strain='Sprage Dawley',
    description='A cute fuzzy animal', # probably don't need this
    sex='M', # The sex of the subject. Using “F” (female), “M” (male), “U” (unknown), or “O” (other) is recommended.
    weight='.3 kg', # Expects kg. Can also be a number instead of string.
)
# If FSCV recordings were performed, use this... In the future, we probably need a custom 'device' and 'data type' for FSCV recordings. For now, we'll just treat it like LFP data.
device = nwbfile.create_device(
    name='FSCV probe',
    description='Carbon probe for measuring FSCV',
    manufacturer='Heien lab'
)

nwbfile.add_electrode_column(name='label', description='label of electrode')

nshanks = 1 # only one carbon probe at a time. 
nchannels_per_shank = 1
electrode_counter = 0

for ishank in range(nshanks):
    # create an electrode group for this shank
    electrode_group = nwbfile.create_electrode_group(
        name='shank{}'.format(ishank),
        description='FSCV probe {}'.format(ishank),
        device=device,
        location='NAc'
    )
    # add electrodes to the electrode table
    for ielec in range(nchannels_per_shank):
        nwbfile.add_electrode(
            x=5.3, y=1.5, z=8.5, imp=np.nan,
            location='NAc',
            filtering='unknown',
            group=electrode_group,
            label='shank{}elec{}'.format(ishank, ielec)
        )
        electrode_counter += 1

nwbfile.electrodes.to_dataframe()

# If LFP recordings were used, use this code...
# - use your data not the random data created here for the demo.
raw_data = np.random.randn(50, 4) # When for real, we will have to put the loaded and processed data here.

all_table_region = nwbfile.create_electrode_table_region(
    region=list(range(electrode_counter)),  # reference row indices 0 to N-1
    description='FSCV all electrodes'
)

# FSCV raw
raw_electrical_series = ElectricalSeries(
    name='ElectricalSeries_FSCV',
    data=raw_data,
    electrodes=all_table_region,
    starting_time=0.,  # timestamp of the first sample in seconds relative to the session start time
    rate=5.  # in Hz
)

nwbfile.add_acquisition(raw_electrical_series)

# LFP - assumes filtered and downsampled. Raw data should be in raw format. I am not sure about how FSCV data should be stored.
lfp_data = np.random.randn(50, 4)
lfp_electrical_series = ElectricalSeries(
    name='ElectricalSeries',
    data=lfp_data,
    electrodes=all_table_region,
    starting_time=0.,
    rate=200.
)
lfp = LFP(electrical_series=lfp_electrical_series)

# Single unit activity. Spike times. This creates some random spike times.
nwbfile.add_unit_column(name='quality', description='sorting quality')

poisson_lambda = 20
firing_rate = 20
n_units = 10
for n_units_per_shank in range(n_units):
    n_spikes = np.random.poisson(lam=poisson_lambda)
    spike_times = np.round(np.cumsum(np.random.exponential(1 / firing_rate, n_spikes)), 5)
    nwbfile.add_unit(spike_times=spike_times,
                     quality='good',
                     waveform_mean=[1., 2., 3., 4., 5.])

nwbfile.units.to_dataframe()

# create some example TimeSeries
test_ts = TimeSeries(name='series1',
                     data=np.arange(1000),
                     unit='m',
                     timestamps=np.linspace(0.5, 601, 1000))
rate_ts = TimeSeries(name='series2',
                     data=np.arange(600),
                     unit='V',
                     starting_time=0.0, rate=1.0)
# Add the TimeSeries to the file
nwbfile.add_acquisition(test_ts)
nwbfile.add_acquisition(rate_ts)

# Intervals 
sleep_stages = TimeIntervals(
    name="sleep_stages",
    description="intervals for each sleep stage as determined by EEG",
)

sleep_stages.add_column(name="stage", description="stage of sleep")
sleep_stages.add_column(name="confidence", description="confidence in stage (0-1)")

sleep_stages.add_row(start_time=0.3, stop_time=0.5, stage=1, confidence=.5)
sleep_stages.add_row(start_time=0.7, stop_time=0.9, stage=2, confidence=.99)
sleep_stages.add_row(start_time=1.3, stop_time=3.0, stage=3, confidence=0.7)

_ = nwbfile.add_time_intervals(sleep_stages) # WTF does the _ do?

nwbfile.get_time_intervals('sleep_stages').to_dataframe()
t = nwbfile.get_acquisition('series1')
plt.figure()
plt.plot(t.data, t.data*0,'.')
plt.xlabel(t.time_unit)
plt.show()
# Write the data: this does work
with NWBHDF5IO('ecephys_tutorial.nwb', 'w') as io:
    io.write(nwbfile)

# Read the data back: this does work
with NWBHDF5IO(r'C:\Users\Stephen Cowen\Documents\GitHub\DANA\src\python\Woody_code\ecephys_tutorial.nwb', 'r') as io:
    read_nwbfile = io.read()
    print(read_nwbfile.acquisition['ElectricalSeries_FSCV'])
    # The following does not work - must be labeling wrong.
    #print(read_nwbfile.processing['ElectricalSeries']) 
    #print(read_nwbfile.processing['ecephys']['LFP'])
    #print(read_nwbfile.processing['ecephys']['LFP']['ElectricalSeries'])