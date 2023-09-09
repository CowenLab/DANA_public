function d = DANA_Load_brain_gel_data_dirs(data_dir, LFP_file,DIN_file,decimation_factor, ...
    template_subtraction, notch_filter_fq,spline_interpolant_removal)

if nargin < 1
    data_dir = pwd;
end
if nargin < 2
    LFP_file = 'amp-D-000.dat';
end
if nargin < 3
    DIN_file = 'board-DIN-03.dat'; % scan file. Probably useless.
end
if nargin < 4
    decimation_factor = 1;
end
if nargin < 5
    template_subtraction = false;
end
if nargin < 6
    notch_filter_fq = []; %make empty if no notch is to be applied, or choose the frequency.
end
if nargin < 7
    spline_interpolant_removal = false;
end
    

d = dir(fullfile(data_dir,'*Hz*'));
d = d([d.isdir]);
for iD = 1:length(d)
    d(iD).name;
    ix = strfind(lower(d(iD).name),'hz');
    if isempty(ix)
        d(iD).Freq = NaN;
        error('This should never happen.')
    else
        d(iD).Freq = str2double(d(iD).name(1:(ix-1)));
    end
end
% Make the very first file the baseline file...
dbase = dir(fullfile(data_dir,'no_input_frequency*'));
dbase(1).Freq = 0;
d(end+1) = dbase(1);
% Re-sort the directories according to ascending frequency.
[~,ix] = sort([d.Freq]);
d = d(ix);
%% Go through each directory and load the data...
for iD = 1:length(d)
    RHD = INTAN_Read_RHD_file(fullfile(data_dir,d(iD).name,'info.rhd'));
    sFreq = RHD.frequency_parameters.board_dig_in_sample_rate;
    DIN_sFreq = sFreq;
    LFP_sFreq = DIN_sFreq/decimation_factor; %
    % Points to ignore from the analysis...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load digital signals.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile(data_dir,d(iD).name,DIN_file),'rb');
    D = fread(fp,'int16');
    fclose(fp);
    %     plot(D(1:1000))
    %     hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the event and block times.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ThreshVal = 0.5;
    EventRecs   =  find(diff(D)>ThreshVal)+1;
    EventRecs_sec = EventRecs/DIN_sFreq;
    EventRecsDown   = find(diff(D)<-ThreshVal)+1;
    if ~isempty(EventRecsDown)
        if EventRecsDown(1) < EventRecs(1);
            EventRecsDown(1) = [];
        end
        if EventRecs(end) > EventRecsDown(end);
            EventRecs(end) = [];
        end
    end
    
    EventRecsDown_sec = EventRecsDown/DIN_sFreq;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the LFP data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile(data_dir,d(iD).name,LFP_file),'rb');
    LFP = fread(fp,'int16');
    fclose(fp);
    LFPunfiltered = LFP;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Denoise- artifact. Useful if the artifact is consistent.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if template_subtraction
        if ~isempty(EventRecs)
            LFP = INTAN_denoise_by_aligned_template_subtraction(LFP,EventRecs(10:2:length(EventRecs)),round(median(diff(EventRecs(1:100)))),EventRecs);
        end
    end
    if spline_interpolant_removal
        if ~isempty(EventRecs)
            pts_before = round(sFreq/50);
            pts_after = round(sFreq/50);
            LFP = INTAN_denoise_by_interpolate_through_artifact_intervals(LFP,[EventRecs(:)-10 EventRecsDown(:)+10],[pts_before pts_after]);
        end
    end
    if ~isempty(notch_filter_fq)
        % This is not a good idea as there are many frequency components to
        % the artifact and notch only gets rid of one very specific
        % frequency.
        [a,b] = INTAN_notch_filter( sFreq, notch_filter_fq, 0.01);
        LFP = filter(b,a,LFP);
    end
    % DO the decimation after all of the template subtraction crap.
    % this reduces a lot of the artifact.
    LFP = decimate(LFP,decimation_factor);
    LFP = single(LFP);
    
    LFP = LFP*RHD.bit_to_uvolt_conversion;

    if length(EventRecs_sec) < length(EventRecsDown_sec)
        EventRecsDown_sec = EventRecsDown_sec(1:length(EventRecs_sec));
    end
    
    x_LFP_sec = (1:length(LFP))'/LFP_sFreq;
    %% Archive in d
    d(iD).LFP = LFP;
    d(iD).x_LFP_sec = x_LFP_sec;
    d(iD).EventRecs_sec = EventRecs_sec;
    d(iD).EventRecsDown_sec = EventRecsDown_sec;
    d(iD).LFP_sFreq = LFP_sFreq;
    d(iD).DIN_sFreq = DIN_sFreq;
    d(iD).RHD = RHD;
    
    %% Plot all of the data and all of the start and end times for a sanity check.
    %
    if false
        figure(1)
        clf
        x_LFPorig_sec = (1:length(LFPunfiltered))'/sFreq;
        plot(x_LFP_sec,LFP,x_LFPorig_sec,LFPunfiltered)
        axis tight
        legend('interp','infiltered')
        hold on
%         plot_markers_simple(EventRecs_sec,[],1,'g')
%         plot_markers_simple(EventRecsDown_sec,[],1,'r')
        ylabel('$\mu$ V','Interpreter','Latex')
        xlabel('Sec')
        pubify_figure_axis
    end
    
end