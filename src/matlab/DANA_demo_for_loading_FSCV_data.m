function OUT = DANA_demo_for_loading_FSCV_data()
%%
% For this to run, you must run when inside a post-processed data
% directory.
% e.g., 
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat425\01
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat445\01
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat439\01
%
% Then just run 
%  DANA_demo_for_loading_FSCV_data
%
% Events.mat
% a folder called WCCV_PCA that has the many .txt tab-separated value files
% and a single .xlsx file that has the order in which the blocks occurred.
% Good sessions:
% GitHub\DANA\Data\Acute\Processed_Data\20221101
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023 (may not work on the 2022 data - change of format since 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
OUT = [];
block_dur_sec = 10;
PLOT_IT = true;
if exist('.\WCCV_PCA','dir')
    WCCV_data_folder = '.\WCCV_PCA';
else
    d = dir('Rat_*_i_vs_t');
    WCCV_data_folder = d(1).name;   
end

%%%%%%%%%%%%%%%%%
% Load the CV data.
[BLOCK_INFO ] = DANA_load_wccv_data_from_folder(WCCV_data_folder);

% To plot data from aa single trial...
% All the meta data for the trial is in BLOCK_INFO
% each element of BLOCK_INFO is a 'trial'
figure
plot(BLOCK_INFO(1).within_block_CV_data(:,1),BLOCK_INFO(1).within_block_CV_data(:,2))
