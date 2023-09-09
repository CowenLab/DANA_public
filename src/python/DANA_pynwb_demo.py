#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NWB Neurodata Without Borders is the format of data for uploading that we need to use.
Once converted, we can store our data (for free) on the DANDI servers.
We can also use a suite of function built for NWB for data analysis should we wish.

Setting things up

Some extra things that you may or may not want
conda install -c conda-forge dandi
pip install nwbwidgets

FOR NEW USERS: 

The only (hopefully) thing you should need to modify are the values in some of the first lines of code that look something like this...
session_data_dir = 'Z:\Data\ACUTE_DANA\220408_MFB_NAc'

nwbfile = NWBFile(
    session_description='FSCV and single unit anesthetized',
    identifier='220408_MFB_NAc',
    .... )

Replace with values appropriate for the experimental session under evaluation.


UPDATE: Currently using Python because the MATLAB version DIES on the first function call.
Also: I found that there are also considerable complications (all due to the stupyd h5py package) with python when you work in environments that have things like spyder installed
I got this to work if I created a fresh CONDA environment (python 3.9) and then installed as follows...
conda install -c conda-forge pynwb. You will also want matplotlib.

Transforming all recording sessions to NWB is required for the Brain Initiative Projects. NWB is foundational for
the Distributed Archives for Neurophysiology Data Integration (DANDI)
data repository to enable collaborative data sharing and analysis.
(https://dandiarchive.org)
Key instructions on data conversion.
https://pynwb.readthedocs.io/en/stable/install_users.html
https://nwb-overview.readthedocs.io/en/latest/conversion_tutorial/

https://www.braininitiative.org/toolmakers/resources/neurodata-without-borders-nwb/
https://www.biorxiv.org/content/10.1101/2021.03.13.435173v2.full

Interfacing: NWBWidgets is a library for interactive web-based visualizations of NWB data, that provides tools for navigating and combining data across modalities.
neurodata types are divided into modules such as ecephys (extracellular
electrophysiology), icephys (intracellular electrophysiology), behavior

an ElectricalSeries is a neurodata type that defines the data and metadata for an intracranially recorded voltage time series in an extracellular electrophysiology experiment.
TimeSeries neurodata type, which is a generic structure designed for any measurement that is sampled over time, and defines fields, such as, data, units of measurement, and sample times (specified either with timestamps or sampling rate and start time)

DANDI can handle large volumes (TBs) of data and host them for free, and
provides a command-line interface (CLI) that is built to handle this
volume of data. DANDI also automatically parses NWB files to extract
metadata that makes these datasets easier for others to find.

https://pynwb.readthedocs.io/en/stable/tutorials/domain/ecephys.htmlsphx-glr-tutorials-domain-ecephys-py


It is easy to set things up wrong. To check, you need to go to the command prompt and run
dandi validate. 
when you are in the directory with the .nwb file.

Some errors I am getting with these data... (I fixed them later - I had multiple channels for the 1 channel FSCV probe by mistake)
2022-07-28 11:47:04,552 [   ERROR]   The second dimension of data does not match the length of electrodes. Your data may be transposed.


To updload to DANDI (UGH)
https://www.dandiarchive.org/handbook/13_upload/

on the command line (not in Python) (be in the directory with the nwb file) 
(get the dataset ID 000299 in this case - will be different for your dataset.)

dandi download https://dandiarchive.org/dandiset/000299/draft
cd 000299
dandi organize .. -f dry
dandi organize ..
dandi upload

It will then ask for the API key. Enter it (from your DANDI site and user icon). A bunch of messages should appear, but it should upload.
When done, if you refresh your DANDI webpage, the set should appear (look in Assets Summary). Also look at Files.

I was able to do all the above and verify that the .nwb file was indeed uploaded.


@author: Cowen

"""
# Required libraries. If it crashes here, it means you did not install the library (using conda or pip)
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
#from nwbwidgets import nwb2widget # useful but not necessary
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
    subject_id = 'Rat203', # PUT IN THE RAT NUMBER (so it can be traced back to the rat order/etc.)
    age='3mo', # for the rat, not the experimenter.
    species='rat', # probably need the propoer scientific name but whatever.
    strain='Sprage Dawley',
    description='A cute fuzzy animal', # probably don't need this
    sex='M', # The sex of the subject. Using “F” (female), “M” (male), “U” (unknown), or “O” (other) is recommended.
    weight='.3 kg', # Expects kg. Can also be a number instead of string.
)
# If FSCV recordings were performed, use this... In the future, we probably need a custom 'device' and 'data type' for FSCV recordings. For now, we'll just treat it like LFP data.
device_FSCV = nwbfile.create_device(
    name='FSCV probe',
    description='Carbon probe for measuring FSCV',
    manufacturer='Heien lab'
)

device_LFP = nwbfile.create_device(
    name='LFP probe',
    description='LFP',
    manufacturer='Cowen Lab'
)

device_Neuropixels = nwbfile.create_device(
    name='Neuropixels probe',
    description='Neuropixels',
    manufacturer='Cowen Lab'
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
        device=device_FSCV, # OK - I think FSCV would be one electrode group. We would create another for the ephys data.
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
raw_data = np.random.randn(50, 1) # When for real, we will have to put the loaded and processed data here.

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
# I am confused - don't I need to also create a 'device' as well for the LFP?
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
    n_spikes = np.random.poisson(lam=poisson_lambda) # Random spiking data.
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