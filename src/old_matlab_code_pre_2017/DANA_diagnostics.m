%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This presumes that we are in the data directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running in \\DISKSTATION\Cowen_Master_Data\Data\DANA\DANA_test_021815\Intan_withSSR_noVolt3_150218_130544
files_to_read = {'amp-D-001.dat' };
%files_to_read = {'amp-D-001.dat' 'amp-D-005.dat' 'amp-D-021.dat'};
analysis_dir = 'C:\Temp';
time_before_event_triggered_average_sec = 0.1;
time_after_event_triggered_average_sec = 0.2;
freq_range_to_plot = [1 500]; % range to plot for the PSD (Pmtm)
event_df_threshold = 1500; % the threshold of the first derivative for event detection.

%
Raw_files = find_files('amp-D*.dat');
IF = INTAN_Read_RHD_file;
Raw_t_sec = INTAN_Load_Time('time.dat'); % 32 bit
sFreq = 1/median(diff(Raw_t_sec));

duration_of_recording_hrs = (Raw_t_sec(end) - Raw_t_sec(1))/(60*60);

% Go through each file of interest
for iF = 1:length(files_to_read)
    D = INTAN_Load_Dat(files_to_read{iF});
    % Detect threshold events
    df = diff(D);
    events_recid = INTAN_Extract_Transitions(files_to_read{iF}, event_df_threshold);
    
    % restrict to only 200 events or less
    if length(events_recid) > 200
        events_recid = events_recid(1:199);
    end
    
    [M, ix] = PETH_EEG_simple([Raw_t_sec(:) double(D(:))], Raw_t_sec(events_recid), ...
        time_before_event_triggered_average_sec*sFreq, time_after_event_triggered_average_sec*sFreq,sFreq);
    x_axis = linspace(-1*time_before_event_triggered_average_sec, time_after_event_triggered_average_sec,Cols(M));
    % Of entire window...
    [pxx,Fqs] = pmtm(double(D(1:40000)),[],[],sFreq);
    % Of a window of time.
    window = 1024; % 64 this seems to be ideal for 1800Hz sFreq ripples.
    overlap = 800; % 56 this seems to be ideal for 1800Hz sFreq ripples.
    [~,fq,T,P] = spectrogram(M(1,:),window,overlap,freq_range_to_plot(1):5:freq_range_to_plot(end),sFreq);
    T = T - time_before_event_triggered_average_sec;
    C = 10*log10(abs(P)); % From the matlab docs.
    [pxx,Fqs] = pmtm(double(D(1:40000)),[],[],sFreq);
    % restrict to the frequencies of interest...
    IX = Fqs > freq_range_to_plot(1) & Fqs < freq_range_to_plot(2);
    pxx = pxx(IX); Fqs = Fqs(IX);
    % do a PSD within one inter-pulse period.
    
    
    % Plot diagnostics.
    figure (1)
    clf
    subplot(3,1,1:2)
    imagesc(T,fq,C)
    xlabel('s')
    ylabel('Hz')
    title(files_to_read{iF})
    subplot(3,1,3)
    plot(x_axis,M(1,:))
    axis tight
    
    [p,n] = fileparts(files_to_read{iF});
    
    saveas(gca,[analysis_dir '\' n '_spectro.png'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(2)
    clf
    subplot(2,2,1:2)
    plot(Raw_t_sec(1:10:end),D(1:10:end))
    xlabel('sec')
    ylabel('? need to figure out how to convert units to ')
    title(files_to_read{iF})

    subplot(2,2,3)
    plot(Fqs,pxx);
    xlabel('Hz')
    ylabel('Power')
    title('Power Spectral Density (PMTM)')
    
    subplot(2,2,4)
    plot_confidence_intervals(x_axis, M);
    title('Event Triggered Average')
    xlabel('sec')
    plot_vert_line_at_zero
    
    saveas(gca,[analysis_dir '\' n '_ETA.png'])

    
end
















