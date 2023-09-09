% This script runs analyses and generates figures used in the 2023 BRAIN poster
% This will people understand how to load and visualze and analyze the data.
%
% Prerequisites: 
% 1) this folder and all subfolders must be in your matlab path
% 2) CowenLib and all subfolders must be in your path.
% - all is on github. Use the most recent version
%
% How I generally organize my code: 
%  general purpose functions are in CowenLib
%  general purpose functions that are still specific to DANA are in the
%  DANA src/matlab folder with the prefix DANA_*.m
%  Analyses that address a specific scientific quesiton are in the
%  ./Questions folder and numbered with Q1_... with Q meaning question and
%  the actual question indicated in the name of the script.
%  Typically, there are 2 types of Questions scripts: the first type Q2*.m
%  should run in the directory that contains the data for a given session.
%  The Q2*_Ana.m aggregates data across sessions and animals to yield an
%  experiment-animal-session-wide analysis.
%
% Cowen 2023
%
% Here are some analyses that went into the BRAIN Initiative poster in
% 2023. They can all run with data stored on GitHub alone (the DANA folder
% and you also need the CowenLib folder) so no need to download a separate
% data file.
%
% The first script (Q2_does_DA_affect_ensemble_act_2023) analyzes datasets 
% where [DA] and neural data were acquired simultaneously and generates 
% some pretty plots that overlay the two signals.
%
% For this to run, you must run when inside a post-processed data
% directory.
% e.g., (the github folder path with differ slightly on your computer)
% Choose one of these...
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat445\01
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat425\01
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat439\01
%
% now that you are in this folder, run the following demos to demontrate
% how to 1) load the FSCV data, and 2) to perform a full analysis that
% integrates both single-unit and dopamine data...

DANA_demo_for_loading_FSCV_data

Q2_does_DA_affect_ensemble_act_2023
% Feel free to copy and rename this function so that you can modify to your
% own needs. I can't say that it is well commented, but it will pop up a
% lot of plots with both ensemble signals and peri-event responses. 

% The second analysis focuses on just the neural ensemble data as there are
% sessions where we have ensemble data but no [DA] data. In addition, this
% script aggregates across multiple recording sessions so analyses a lot of
% data. You do not need to be in a specific subfolder for this function to
% work (unlike the one above). The function below relies on data processed
% in each session folder by Q3_does_stim_affect_ensemble_act.m. I already
% ran this and stored the output on GitHub, but if you would like to know
% how the data used by Q3_does_stim_affect_ensemble_act_Ana.m were
% generated, look at Q3_does_stim_affect_ensemble_act.m
data_dir = 'C:\Users\cowen\Documents\GitHub\DANA\Data\Q3_does_stim_affect_ensemble_act';
% Again, the above will change slightly depending on the path on your
% computer to the DANA Github folder.
% Now let's run the script.
Q3_does_stim_affect_ensemble_act_Ana(data_dir)
% Feel free to copy and rename this function so that you can modify to your
% own needs.




