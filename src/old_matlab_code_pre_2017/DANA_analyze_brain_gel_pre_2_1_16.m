function [OUT,META] = DANA_analyze_brain_gel_pre_2_1_16(data_dir)
% Create brain gel figures for the DANA project.
%  - assumes INTAN support functions are in the path.
% Assumes that you are in a directory with the rhd, dat, and datlfp data.
%
% NOTE: This works on pre 2/24/2016 files. Things changed and a new version
% was created to handle more recent unlabeled data.
% 
% Cowen 2015 (with some inspiration from JP's code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load meta data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = []; META = [];
if nargin < 1
    data_dir = pwd;
end
%%
mname = 'DANA_analyze_brain_gel';
RHD = INTAN_Read_RHD_file(fullfile(data_dir,'info.rhd'));
DIN_sFreq = RHD.frequency_parameters.board_dig_in_sample_rate;
LFP_sFreq = DIN_sFreq/2; %
% Points to ignore from the analysis...
META.Time_from_start_to_ignore_sec = 0.01;
pts_to_ignore = META.Time_from_start_to_ignore_sec *LFP_sFreq; % ignore 5 
sec_before = 10; % time before to determine spectral power.
Freq_high =  RHD.frequency_parameters.actual_upper_bandwidth;
Attenuation_table = [.5*Freq_high .99; .8*Freq_high  0.89; Freq_high  0.707]; % from the datasheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIN_file = 'board-DIN-03.dat'; % scan file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP_file = 'LFP_Ch_3_Skip_10.datlfp';skip_factor = 10;
LFP_file = 'amp-D-000.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load digital signals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile(data_dir,DIN_file),'rb');
D = fread(fp,'int16');
fclose(fp);

if 0
    d = dir('board-DIN*.dat');
    for ii =1:length(d)
        fp = fopen(fullfile(data_dir,d(ii).name),'rb');
        D = fread(fp,'int16');
        fclose(fp);
        d(ii).name
        sum(D>0)
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the event and block times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThreshVal = 0.5;
EventRecs   =  find(diff(D)>ThreshVal)+1;
EventRecs_sec = EventRecs/DIN_sFreq;
EventRecsDown   = find(diff(D)<-ThreshVal)+1;
EventRecsDown_sec = EventRecsDown/DIN_sFreq;
% Block Times
BlockEnds_sec = [EventRecs_sec(diff(EventRecs_sec)>20); EventRecs_sec(end)];
BlockStarts_sec = EventRecs_sec(diff([0;EventRecs_sec])>20);
Block_frequencies = [0 1 2 5 7 10 14 20 40 50 80 100 140 200 500 1000 2000 2857];
Frequency_range_to_examine = [0 0 40;
    1 0 40;
    2 0 40;
    5 0 40;
    7 0 40;
    10 0 80;
    14 0 80;
    20 0 100;
    40 0 100;
    50 0 100;
    80 10 150;
    100 10 150;
    140 20 200;
    200 40 300;
    500 80 700;
    1000 600 1800;
    2000 1500 2400;
    2857 2400 3600; ];

Frequency_for_SigToNoise = [0 0 0;
    1 .5 5;
    2 .5 6;
    5 1 10;
    7 3 15;
    10 5 20;
    14 8 22;
    20 10 30;
    40 25 55;
    50 35 60;
    80 65 95;
    100 80 120;
    140 120 160;
    200 180 240;
    500 400 600;
    1000 800 1200;
    2000 1800 2200;
    2857 2600 3100; ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the LFP data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile(data_dir,LFP_file),'rb');
LFP = fread(fp,'int16');
fclose(fp);
LFP = LFP*RHD.bit_to_uvolt_conversion;
%% Plot all of the data and all of the start and end times for a sanity check.
x_LFP_sec = (1:length(LFP))'/LFP_sFreq;
%
if 0
    figure(1)
    clf
    plot(x_LFP_sec/60,LFP)
    axis tight
    hold on
    plot_markers_simple(EventRecs_sec/60,[],[],'g')
    plot_markers_simple(EventRecsDown_sec/60,[],[],'r')
    plot_markers_simple(BlockStarts_sec/60,[],[],'b')
    plot_markers_simple(BlockEnds_sec/60,[],[],'m')
end
% Get response during baseline...
LFP_baseline = LFP(x_LFP_sec > (BlockStarts_sec(1) - sec_before) & x_LFP_sec < (BlockStarts_sec(1) - 2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each block, generate the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iBlock = 2:length(BlockStarts_sec)
    ix_before = find(x_LFP_sec > (BlockStarts_sec(iBlock) - sec_before) & x_LFP_sec < (BlockStarts_sec(iBlock) - 2) );
    ix_after = find(x_LFP_sec > BlockStarts_sec(iBlock) & x_LFP_sec < BlockEnds_sec(iBlock)  );
    EventRecs_sec_block = EventRecs_sec(EventRecs_sec>=BlockStarts_sec(iBlock) & EventRecs_sec<=BlockEnds_sec(iBlock));
    EventRecsEnds_sec_block = EventRecsDown_sec(EventRecsDown_sec>EventRecs_sec_block(1) & EventRecsDown_sec<=BlockEnds_sec(iBlock));
    if EventRecs_sec_block(end) > EventRecsEnds_sec_block(end)
        EventRecs_sec_block(end) = [];
    end
    EventRecs_sec_block = EventRecs_sec_block - BlockStarts_sec(iBlock);
    EventRecsEnds_sec_block  = EventRecsEnds_sec_block - BlockStarts_sec(iBlock);
    welch_fq = Frequency_range_to_examine(iBlock,2):1:Frequency_range_to_examine(iBlock,3);
    SnR_FQ = [Frequency_for_SigToNoise(iBlock,2) Block_frequencies(iBlock) Frequency_for_SigToNoise(iBlock,3)];

    %     time_to_plot_sec = [0.005 median(diff(EventRecs_sec_block))-.005];
    
    LFP_before = LFP(ix_before);
    LFP_after = LFP(ix_after);
    %     x_sec_before = x_LFP_sec(ix_before)-x_LFP_sec(ix_before(1));
    x_sec_after = x_LFP_sec(ix_after)-x_LFP_sec(ix_after(1));
    
    if 0
        figure
        plot(x_sec_after(:),LFP_after(:))
        hold on
        plot_markers_simple(EventRecs_sec_block,[],[],'g')
        plot_markers_simple(EventRecsEnds_sec_block)
    end
    window_sec = min([20*1/Block_frequencies(iBlock) 5 length(LFP_before)/LFP_sFreq]);
    window_samples_spectrogram = round(window_sec * LFP_sFreq);
    %     window_samples_welch = 10+round((LFP_sFreq*2.5)./Block_frequencies(iBlock));
    
    %     window_samples = 2*2^14;
    
    [PW_baseline] = pwelch(LFP_baseline,[],[],welch_fq,LFP_sFreq); % let welch decide.
    [PW_snr_baseline] = pwelch(LFP_baseline,[],[],SnR_FQ,LFP_sFreq); % let welch decide.
    
    [~,fq,~,P] = spectrogram(LFP_before,window_samples_spectrogram,round(window_samples_spectrogram/2),[],LFP_sFreq);
    S_before = 10*log10(abs(P)); % From the matlab docs.
    x_S_sec = linspace(x_LFP_sec(ix_before(1)),x_LFP_sec(ix_before(end)),Cols(S_before));
    fq_IX = fq > Frequency_range_to_examine(iBlock,2) & fq <Frequency_range_to_examine(iBlock,3);
    S_before = S_before(fq_IX,:);
    S_fq = fq(fq_IX);
    [PW_before] = pwelch(LFP_before,[],[],welch_fq,LFP_sFreq); % let welch decide.
    [PW_snr_before] = pwelch(LFP_before,[],[],SnR_FQ,LFP_sFreq); % let welch decide.

    %[PW, welch_fq] = pwelch(LFP_before,window_samples_welch,round(window_samples_welch/2),[],LFP_sFreq);
    %[PW, welch_fq] = pmtm(LFP_before,window_samples,[],LFP_sFreq); DO NOT
    %USE PMTM - at least how I implemented it. Too slow.
    % Units PSD estimate, returned as a real-valued, nonnegative column vector or matrix. Each column of pxx is the PSD estimate of the corresponding column of x. The units of the PSD estimate are in squared magnitude units of the time series data per unit frequency. For example, if the input data is in volts, the PSD estimate is in units of squared volts per unit frequency. For a time series in volts, if you assume a resistance of 1 ? and specify the sampling frequency in hertz, the PSD estimate is in watts per hertz.
    pwelch_before = 10*log10(PW_before);
    [~,ix] = max(pwelch_before);
    pq_fq = welch_fq(ix);
    
    ixfqlow = find(welch_fq >= Frequency_for_SigToNoise(iBlock,2),1,'first');
    ixfqhigh = find(welch_fq >= Frequency_for_SigToNoise(iBlock,3),1,'first');
    ixfqtgt = find(welch_fq >= Frequency_for_SigToNoise(iBlock,1),1,'first');
    SnR_before_welch = PW_before(ixfqtgt)/mean(PW_before([ixfqlow ixfqhigh]));
    
    OUT(iBlock).Target_Fq = Block_frequencies(iBlock);
    OUT(iBlock).SnR_before_welch = SnR_before_welch;
    OUT(iBlock).welch_fq = welch_fq;
    OUT(iBlock).pwelch_before = pwelch_before;
    OUT(iBlock).Frequency_for_SigToNoise = Frequency_for_SigToNoise(iBlock,:);
    OUT(iBlock).Peak_mV_from_hilbert_before = mean(abs(hilbert(LFP_before)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the frequency response during baseline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(iBlock*10)
    clf
    subplot(3,1,1)
    plot(x_LFP_sec(ix_before) - x_LFP_sec(ix_before(1)), LFP(ix_before))
    axis tight
    plot_horiz_line_at_zero(OUT(iBlock).Peak_mV_from_hilbert_before)

    xlabel('s'); ylabel('uV')
    pubify_figure_axis
    title([num2str(Block_frequencies(iBlock)) ' Hz Baseline'])
    
    subplot(3,1,2)
    cla
    imagesc(x_S_sec-x_S_sec(1),S_fq,S_before)
    %     set(gca,'XTick',0:(Cols(S_before)-1))
    colormap(jet)
    axis xy
    set(gca,'XTickLabel',[])
    xlabel('s'); ylabel('Frequency (Hz)')
    title('Spectrogram')
    pubify_figure_axis
    colorbar_label('log power')
    
    subplot(3,1,3)
    %plot(S_fq, mean(S_before,2))
    plot(welch_fq, pwelch_before)
    xlabel('Frequency (Hz)');ylabel('10*log10(squared uV/Hz)')
    axis tight
    title(sprintf('Welch PSD, Max Fq = %2.4f, SnR = %2.4f SnR Border Fq %2.1f %2.1f Hz ',pq_fq, SnR_before_welch, Frequency_for_SigToNoise(iBlock,2), Frequency_for_SigToNoise(iBlock,3)))
    pubify_figure_axis
    label_figure(mname);
    % Save figure.
    saveas(gcf,fullfile('Figures',sprintf('Before_Fq_%f.png',Block_frequencies(iBlock))),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the ratio of the power in the outer frequencies compared to the
    % target frequency for each band.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First, iterate through each interval and determine the spectrogram
    % for each inter-pulse period.
    dur_s = median(EventRecs_sec_block(2:end)-EventRecsEnds_sec_block(1:(end-1)));
    newx = linspace(0, dur_s,round(dur_s*LFP_sFreq));
    welch_IX = newx > META.Time_from_start_to_ignore_sec; % & newx < newx(end)-0.01;
    spec_IX = newx > META.Time_from_start_to_ignore_sec; % & newx < newx(end)-0.01;
    L = zeros(length(newx), length(EventRecs_sec_block)-1);
    last_ix = round(length(newx)-(0.05*LFP_sFreq):length(newx));
    for iScan = 1:(length(EventRecs_sec_block)-1)
        IX = x_sec_after >= EventRecsEnds_sec_block(iScan) & x_sec_after <= EventRecs_sec_block(iScan+1);
        x = x_LFP_sec(IX) - EventRecs_sec_block(iScan);
        L(:,iScan) = interp1(x,LFP_after(IX),newx);
        L(~spec_IX,:) = nan;
        NANIX = ~isnan(L(:,iScan));
        L(NANIX,iScan) = detrend(L(NANIX,iScan));
        L(:,iScan) = L(:,iScan) - nanmean(L(last_ix,iScan));
    end
    
    Frequency_for_SigToNoise(iBlock,2:3)
    b = fir1(95,Frequency_for_SigToNoise(iBlock,2:3)*2/LFP_sFreq);
%     [aa, bb ]= butter(5,Frequency_for_SigToNoise(iBlock,2:3)*2/LFP_sFreq);
 
    if 0
        b = fir1(148,Frequency_for_SigToNoise(iBlock,2:3)*2/LFP_sFreq);
        [aa, bb ]= butter(5,Frequency_for_SigToNoise(iBlock,2:3)*2/LFP_sFreq);
        figure
        freqz(b,1,512,LFP_sFreq)
        figure
        freqz(bb,aa,512,LFP_sFreq)
        LFP_low = filter(b,1,LFP(ix_before));
    end
    % We need good time resolution...
    %     window_sec = min([5/Block_frequencies(iBlock) 0.05]); %  a miniumum of around 5 cycles seems like a bare minium
    window_samples = min([3*LFP_sFreq./Block_frequencies(iBlock) 0.02*LFP_sFreq]);
    window_samples = round(window_samples);
    [~,fq,~,P] = spectrogram(LFP_before,window_samples,round(window_samples/3),SnR_FQ,LFP_sFreq);
    S_before = abs(P); %
    %     S_before = 10*log10(abs(P)); % From the matlab docs.
    x_S_before_sec = linspace(x_LFP_sec(ix_before(1)),x_LFP_sec(ix_before(end)),Cols(S_before));
    SnR_before = S_before(2,:)./mean(S_before([1 3],:));
    SNR_by_time = [];
    Peak_mV_from_hilbert_after = zeros(Cols(L),1);
    Peak_mV_from_hilbert_after2 = zeros(Cols(L),1);
    Lv = L.*0;
    Lvf = L.*0;
    
    for iScan = 1:Cols(L)
        v = detrend(L(welch_IX,iScan));
        v_filt1 = filter(b,1,v);
        %v_filt1 = filfilt(bb,aa,v);
        mid = round(length(v)/2);
        Lv(welch_IX,iScan) = v;
        Lvf(welch_IX,iScan) = v_filt1;
        
        [~ , fq, ~, P]  = spectrogram( v, window_samples,...
            round( window_samples/3), SnR_FQ, LFP_sFreq);
        if iScan == 1
            [~ , tfq, ~, tmp]  = spectrogram( v, window_samples,...
                round( window_samples/3), welch_fq, LFP_sFreq);
            fullP = zeros(size(tmp,1),size(tmp,2),Cols(L))*nan;
            fullP(:,:,iScan) = tmp;
        else
            [~ , tfq, ~, fullP(:,:,iScan)]  = spectrogram( v, window_samples,...
                round( window_samples/3), welch_fq, LFP_sFreq);
        end
        
        [PW] = pwelch( v, sum(welch_IX), [], welch_fq, LFP_sFreq);
        PWn = PW - PW_baseline;
        [PWsnr, ~] = pwelch(v, sum(welch_IX), [],SnR_FQ , LFP_sFreq);
        PWsnr_n = PWsnr - PW_snr_baseline;

        if iScan == 1
            S_after_all = zeros(Rows(P),Cols(P), Cols(L));
            pwelch_after = zeros(length(PW),Cols(L));
            pwelch_after_snr = zeros(length(PWsnr),Cols(L));
            pwelch_after_norm = zeros(length(PW),Cols(L));
            pwelch_after_snr_norm = zeros(length(PWsnr),Cols(L));
        end
        S_after_all(:,:,iScan) = abs(P); % From the matlab docs.
        %         S_after_all(:,:,iScan) = 10*log10(abs(P)); % From the matlab docs.
        % Units PSD estimate, returned as a real-valued, nonnegative column vector or matrix. Each column of pxx is the PSD estimate of the corresponding column of x. The units of the PSD estimate are in squared magnitude units of the time series data per unit frequency. For example, if the input data is in volts, the PSD estimate is in units of squared volts per unit frequency. For a time series in volts, if you assume a resistance of 1 ? and specify the sampling frequency in hertz, the PSD estimate is in watts per hertz.
        pwelch_after(:,iScan) = abs(PW)';
        pwelch_after_snr(:,iScan) = abs(PWsnr)';
        pwelch_after_norm(:,iScan) = abs(PWn)';
        pwelch_after_snr_norm(:,iScan) = abs(PWsnr_n)';
        %         pwelch_after(:,iScan) = 10*log10(PW)';
        %         pwelch_after_snr(:,iScan) = 10*log10(PWsnr)';
        % Determine how long it takes for the SNR of 1.5 to be reached...
        range_to_test = 100:50:length(v);
        PWsnrT = zeros(length(range_to_test),3)*nan;
        for ii = 1:length(range_to_test)
            ix = 1:range_to_test(ii);
            PWsnrT(ii,:) = pwelch(v(ix), length(ix), [],SnR_FQ , LFP_sFreq);
        end
        SNR_by_time(:,iScan) = PWsnrT(:,2)./mean(PWsnrT(:,[1 3]),2);
        Peak_mV_from_hilbert_after(iScan) = mean(abs(hilbert(v(mid:end))));
        Peak_mV_from_hilbert_after2(iScan) = mean(abs(hilbert(v_filt1(mid:end))));
    end
    figure
    clf
    plot(mean(L'))
    hold on
    plot(mean(Lv'))
    plot(mean(Lvf'))
    legend('L','Lv','Lvf')
    
    
    if 0 
        plot_confidence_intervals(newx(range_to_test),SNR_by_time')

    end
    S_after = nanmean(S_after_all,3);
    
    x_S_after_sec = linspace(newx(1),newx(end),Cols(S_after));
    SnR_after = S_after(2,:)./mean(S_after([1 3],:));
    SnR_after_welch = pwelch_after_snr(2,:)./mean(pwelch_after_snr([1 3],:));
    SnR_after_welch_norm = pwelch_after_snr_norm(2,:)./mean(pwelch_after_snr_norm([1 3],:));
    
    [~,ix] = max(mean(pwelch_after,2));
    pq_fq = welch_fq(ix);
    
    OUT(iBlock).Mean_Welch_Power_after = mean(pwelch_after_snr,2)';
    OUT(iBlock).Mean_Welch_Power_after_norm = mean(pwelch_after_snr_norm,2)';
    OUT(iBlock).Mean_Welch_SnR_after = mean(SnR_after_welch);
    OUT(iBlock).Mean_Welch_SnR_after_norm = mean(SnR_after_welch_norm);
    OUT(iBlock).Peak_Fq_after = pq_fq;
    OUT(iBlock).Mean_Full_Welch_after = mean(pwelch_after,2);
    OUT(iBlock).SNR_by_time =mean(SNR_by_time,2);
    OUT(iBlock).SNR_by_time_x_sec =newx(range_to_test);
    OUT(iBlock).Peak_mV_from_hilbert_after = median(Peak_mV_from_hilbert_after);
    %%%%%%%%%%%%%%%%%%%%
    figure(iBlock*100)
    clf
    subplot(3,2,1)
    plot(newx,L)
    axis tight; xlabel('s');ylabel('uV'); pubify_figure_axis
    hold on
    plot(0,OUT(iBlock).Peak_mV_from_hilbert_after,'r>');
    a = axis; a(1:2) = [0 dur_s]; axis(a);
    title([num2str(Block_frequencies(iBlock)) ' Hz'])
    
    pubify_figure_axis
    subplot(3,2,2)
    
    plot(x_S_after_sec,S_after,'o-')
    axis tight; xlabel('s');ylabel('Power'); pubify_figure_axis
    title(sprintf('Max Fq = %f Hz',pq_fq))
    legend('lower','target','upper')
    
    subplot(3,2,3)
    % Spectrogram.
    cla
    x = linspace(x_S_after_sec(1) + META.Time_from_start_to_ignore_sec,x_S_after_sec(end),size(fullP,2));
    imagesc(x, tfq,mean(fullP,3))
    colormap('jet')
    axis tight
    axis xy
    ylabel('Hz')
    title('Power')
    colorbar
    %   cla
    %     plot(x_S_before_sec - x_S_before_sec(1), SnR_before)
    %     axis tight; xlabel('s');ylabel('SnR'); pubify_figure_axis
    %     title('SnR Before....?')
    %     plot_horiz_line_at_zero(1)
    
    subplot(3,2,4)
    plot(x_S_after_sec, SnR_after,'o-')
    axis tight; xlabel('s');ylabel('SnR'); pubify_figure_axis
    title('SnR After Spectrogram')
    plot_horiz_line_at_zero(1)
    
    
    subplot(3,2,5)
    plot_confidence_intervals(welch_fq, pwelch_after')
    
    axis tight; xlabel('Freq (Hz)');ylabel('Welch PSD Power'); pubify_figure_axis
    
    subplot(3,2,6)
    cla
    errorb(1:2, nanmean([SnR_after_welch(:) SnR_after_welch_norm(:) ]), Sem([SnR_after_welch(:) SnR_after_welch_norm(:)]))
    hold on
    
    plot_horiz_line_at_zero(1);
    ylabel('SnR'); pubify_figure_axis
    set(gca,'XTickLabel',{'SnR After' 'SnR After Norm'})
    title('From total Welch PSD')
    label_figure(mname);
    
    % Save figure.
    saveas(gcf,fullfile('Figures',sprintf('Fq %f_summary.png',Block_frequencies(iBlock))),'png')
    
    
    fprintf('.')
end
% Save pre data.
save(mname);
%%
SnR = [OUT.Mean_Welch_SnR_after];
FQ  = [OUT.Target_Fq];
xIX = FQ < 100;

figure
bar((SnR(xIX)),'k')
set(gca,'XTick',1:length(SnR(xIX)))
set(gca,'XTickLabel',FQ(xIX))
plot_horiz_line_at_zero();
pubify_figure_axis
ylabel('SnR')
xlabel('Frequency (Hz)')
grid on
label_figure(mname);

