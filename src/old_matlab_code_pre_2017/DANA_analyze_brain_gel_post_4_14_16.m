function [OUT,META] = DANA_analyze_brain_gel_post_4_14_16()
% Create brain gel figures for the DANA project.
%  - assumes INTAN support functions are in the path.
% Assumes that you are in a directory with the rhd, dat, and datlfp data.
% or pass in the data_dir with the .dat files:
% e.g., data_dir = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF on';
%
%
% Cowen 2015 (with some inspiration from JP's code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load meta data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'defaultAxesFontName', 'Arial');
mname = 'DANA_analyze_brain_gel_post_4_14_16';
saveit = false;
without_SSR = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP_file = 'LFP_Ch_3_Skip_10.datlfp';skip_factor = 10;
if without_SSR
    data_dir_wfon = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF on';
    data_dir_wfoff = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF off';
    condition = 'No SSR';
else
    data_dir_wfon = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\with SSR WF on';
    data_dir_wfoff = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\with SSR  WF off';
    condition = 'With SSR';
end
data_dir_wfoff_no_ssr = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF off';
% This is used to generatee our idealized values of the power spectra.

OUT.condition = condition;
subtract_template = false; % this works for high-frequency stuff - not so great for low frequency as teh template amplitudeis not smart enough to adapt to the baseline.
interpolate_out_artifact = false;
% do not do this right now. Not so great with slow frequency artifact. Nice
% for HF noise though. Maybe post-filtering.
%%
sec_before = 30;
wavelet_size = 20;
wavelet_size_in_scan = wavelet_size; % These two need to be the same for now. Problem - values probably need to change radically depending on the frequency - lower for low frequencies, higher for high frequencies.

window_sec =.01;
welch_win = [];
FOR_HIGH_FREQUENCIES = false;
META.Time_from_start_to_ignore_sec = 0.01;
META.Time_from_start_to_ignore_sec = 0.002;
META.fq_logspace = unique(ceil(logspace(log10(1),log10(3000),200)));
% d = [1 diff(META.fq_logspace)];
% for ii = 2:length(d)
% end
% META.fq_logspace = META.fq_logspace(d > 2); % 2 Hz differences are uninteresting.


decimation_factor = 3; % set to 1 for NO decimation
OUT.META = META;
txt = sprintf('%s_tmpsub%d_intrpart%d',condition,subtract_template,interpolate_out_artifact);
mname = [mname txt];

PLOT_IT = true;
OUT = [];
if nargin < 1
    data_dir = pwd;
end
if ~exist(fullfile(data_dir,'Figures'),'dir')
    mkdir(fullfile(data_dir,'Figures'))
end
%% Create a list of sub directories.
don = DANA_Load_brain_gel_data_dirs(data_dir_wfon, 'amp-D-000.dat', 'board-DIN-03.dat',decimation_factor,subtract_template,[],interpolate_out_artifact);
doff = DANA_Load_brain_gel_data_dirs(data_dir_wfoff, 'amp-D-000.dat', 'board-DIN-03.dat',decimation_factor,false);
doff_no_ssr = DANA_Load_brain_gel_data_dirs(data_dir_wfoff_no_ssr, 'amp-D-000.dat', 'board-DIN-03.dat',decimation_factor,false);
%%
% Block Times
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
    1 .2 7;
    2 .5 9;
    5 1 9;
    7 3 11;
    10 7 17;
    14 7 23;
    20 13 37;
    40 33 47;
    50 43 67;
    80 67 97;
    100 83 117;
    140 123 163;
    200 183 243;
    500 403 603;
    1000 803 1203;
    2000 1803 2203;
    2857 2603 3103; ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get response during baseline...
% base_st_ed_sec = Events.FqEpochs(1,2:3);
IX_BASELINE = (don(1).x_LFP_sec > don(1).x_LFP_sec(end) - sec_before) & (don(1).x_LFP_sec < don(1).x_LFP_sec(end) - 4);
LFP_baseline = don(1).LFP(IX_BASELINE);
% LFPsub_baseline = don(1).LFPsub(IX_BASELINE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each block, generate the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iD = 1:length(don)
    tit = sprintf('%s %d Hz',condition, don(iD).Freq);
    block_fq = don(iD).Freq;
    LFP_sFreq = don(iD).LFP_sFreq;
    %     window_samples = max([1.5*LFP_sFreq./block_fq 0.02*LFP_sFreq]);
    if block_fq == 0
        window_samples = (1.5*LFP_sFreq./1);
    else
        window_samples = (1.5*LFP_sFreq./block_fq);
        window_samples = ceil(window_samples);
    end
    ix = find(Frequency_for_SigToNoise(:,1) == block_fq);
    Frequency_for_SigToNoise_BLK = Frequency_for_SigToNoise(ix,:);
    ix = find(Frequency_range_to_examine(:,1) == block_fq);
    Frequency_range_to_examine_BLK = Frequency_range_to_examine(ix,:);
    welch_fq = Frequency_range_to_examine_BLK(2):.5:Frequency_range_to_examine_BLK(3);
%     S_fq = Frequency_range_to_examine_BLK(2): Frequency_range_to_examine_BLK(3);
    S_fq = welch_fq;

    SnR_FQ = [Frequency_for_SigToNoise_BLK(2) block_fq Frequency_for_SigToNoise_BLK(3)];
    if PLOT_IT
        IX_ON = don(iD).x_LFP_sec >= 5 & don(iD).x_LFP_sec <= 5.5;
        figure(1)
        clf
        plot(don(iD).x_LFP_sec(IX_ON), don(iD).LFP(IX_ON))
        hold on
        plot(doff(iD).x_LFP_sec(IX_ON), doff(iD).LFP(IX_ON))
        legend('WF On','WF Off'); legend boxoff
        title(tit)
        pubify_figure_axis
        ylabel('$\mu$V','Interpreter','Latex')
        xlabel('sec')
        if saveit, saveas(gcf,sprintf('SinglesubplotExample%s.eps',tit)), end

    end
    %    LFP_off = LFP(ix_before);
    %    LFP_on  = LFP(ix_after);
    %     x_sec_before = x_LFP_sec(ix_before)-x_LFP_sec(ix_before(1));
    IXON = don(iD).x_LFP_sec >= don(iD).x_LFP_sec(end) - 20 & don(iD).x_LFP_sec <= don(iD).x_LFP_sec(end) - 5;
    IXOFF = doff(iD).x_LFP_sec >= doff(iD).x_LFP_sec(end) - 20 & doff(iD).x_LFP_sec <= doff(iD).x_LFP_sec(end) - 5;
    IXOFF_nossr = doff_no_ssr(iD).x_LFP_sec >= doff_no_ssr(iD).x_LFP_sec(end) - 20 & doff_no_ssr(iD).x_LFP_sec <= doff_no_ssr(iD).x_LFP_sec(end) - 5;
    
    if PLOT_IT
        IXON2 = don(iD).x_LFP_sec >= don(iD).x_LFP_sec(1) + 10 & ...
            don(iD).x_LFP_sec <= 11;
        IXOFF2 = doff(iD).x_LFP_sec >= doff(iD).x_LFP_sec(1) + 10 & ...
            doff(iD).x_LFP_sec <= 11;

        figure
%         L1 = 
        
        DANA_plot_wavelet_and_scans( don(iD).LFP(IXON2),don(iD).x_LFP_sec(IXON2),...
            LFP_sFreq,don(iD).EventRecsDown_sec,S_fq,wavelet_size_in_scan);
        hold on
        plot(doff(iD).x_LFP_sec(IXOFF2),doff(iD).LFP(IXOFF2)/1000,'c')
        title(tit)
        
            figure
%         L1 = 
        
        DANA_plot_wavelet_and_scans( don(iD).LFP(IXON2),don(iD).x_LFP_sec(IXON2),...
            LFP_sFreq,don(iD).EventRecsDown_sec,META.fq_logspace ,wavelet_size_in_scan);
        hold on
        plot(doff(iD).x_LFP_sec(IXOFF2),doff(iD).LFP(IXOFF2)/1000,'c')
        title(tit)
        
    end
    %     window_sec = min([5 sum(IXON)/LFP_sFreq]);
    
    %     window_samples_spectrogram = round(window_sec * LFP_sFreq);
    %     window_samples
    %     window_samples_welch = 10+round((LFP_sFreq*2.5)./block_fq);
    
    %     window_samples = 2*2^14;
    if subtract_template
        LFP_on = don(iD).LFPsub(IXON);
    else
        LFP_on = don(iD).LFP(IXON);
    end
    OUT(iD).subtract_template = subtract_template;
    
    LFP_on_x_sec =don(iD).x_LFP_sec(IXON);
    LFP_on_x_sec = LFP_on_x_sec - LFP_on_x_sec(1);
    
    sFreq = LFP_sFreq;
    
    LFP_off = doff(iD).LFP(IXOFF);
    LFP_off_x_sec =doff(iD).x_LFP_sec(IXOFF);
    LFP_off_x_sec = LFP_off_x_sec - LFP_off_x_sec(1);
        
    LFP_off_no_ssr = doff_no_ssr(iD).LFP(IXOFF_nossr);
    LFP_off_no_ssr_x_sec = doff_no_ssr(iD).x_LFP_sec(IXOFF_nossr);
    LFP_off_no_ssr_x_sec = LFP_off_no_ssr_x_sec - LFP_off_no_ssr_x_sec(1);

    
    [PW_FSCVOFF] = pwelch(LFP_off,welch_win,[],welch_fq,sFreq); % let welch decide.

    [PW_FSCVOFF_noSSR] = pwelch(LFP_off_no_ssr,welch_win,[],welch_fq,sFreq); % let welch decide.
    
    PW_before = PW_FSCVOFF; % remember that this is already log scaled.
    PW_before_no_SSR = PW_FSCVOFF_noSSR; % remember that this is already log scaled.

    [PW_snr_before] = pwelch(LFP_off,welch_win,[],SnR_FQ,sFreq); % let welch decide.
    SnR_before_welch = PW_snr_before(2)/mean(PW_snr_before([1 3]));
    
    [PW_FSCVON] = pwelch(LFP_on,welch_win,[],welch_fq,sFreq); % let welch decide.
    PW_after = PW_FSCVON; % remember that this is already log scaled.
    [SnR_PWR_after_welch_all] = pwelch(LFP_on,welch_win,[],SnR_FQ,sFreq); % let welch decide.
    SnR_after_welch_all = SnR_PWR_after_welch_all(2)/mean(SnR_PWR_after_welch_all([1 3]));
    SnR_r_val_after = corrcoef(PW_before_no_SSR(:),PW_after(:));
    SnR_r_val_after = SnR_r_val_after(2);

    
    if PLOT_IT
        IX = LFP_off_x_sec< 1;
        figure(4)
        clf
        subplot(3,1,1)
        plot(LFP_off_x_sec(IX), LFP_off(IX))
        title([tit ' FSCV Off'])
        ylabel('$\mu$V','Interpreter','Latex')
        
        axis tight
        xlabel('s')
        subplot(3,1,2)
        plot(LFP_on_x_sec(IX), LFP_on(IX))
        axis tight
        title('FSCV On')
        xlabel('s')
        ylabel('$\mu$V','Interpreter','Latex')
        
        subplot(3,1,3)
        plot(welch_fq, 10*log10(PW_FSCVOFF),'LineWidth',2)
        hold on
        plot(welch_fq, 10*log10(PW_FSCVON),'LineWidth',2)
        plot(welch_fq, 10*log10(PW_FSCVOFF_noSSR),'LineWidth',2)

        legend('FSCVOFF SIGOFF','FSCVON SIGOFF','OFFNoSSR')
        xlabel('Hz'); axis tight
        ylabel('Power dB/Hz (10log10)','Interpreter','Latex')
        %         pubify_figure_axis
        if saveit, saveas(gcf,sprintf('Example%s.eps',tit)), end

    end
    
    [~,fq,~,P] = spectrogram(LFP_off,window_samples,round(window_samples/2),S_fq,sFreq);
    S_before = 10*log10(abs(P)); % From the matlab docs. imagesc(fq,[],S_before)
    x_S_sec = linspace(LFP_off_x_sec(1),LFP_off_x_sec(end),Cols(S_before));

    
    wavelet_fqs = fq;
    [~,Pwv] = SPEC_waveletdecomp(S_fq,LFP_off,sFreq,wavelet_size);
    Pwv = Pwv(:,round(2*sFreq):(end-round(2*sFreq)));
    
    [~,Pwv_no_ssr] = SPEC_waveletdecomp(S_fq,LFP_off_no_ssr,sFreq,wavelet_size);
    Pwv_no_ssr = Pwv_no_ssr(:,round(2*sFreq):(end-round(2*sFreq)));
    
    [~,Pwv_on] = SPEC_waveletdecomp(S_fq,LFP_on,sFreq,wavelet_size);
    Pwv_on = Pwv_on(:,round(2*sFreq):(end-round(2*sFreq)));
    


    MIDIX = round(size(Pwv_no_ssr,2)/2-100:size(Pwv_no_ssr,2)/2+100);
    [~,ix] = max(PW_after);
    pq_fq = welch_fq(ix);

    
    OUT(iD).Target_Fq = block_fq;
    OUT(iD).SnR_before_welch = SnR_before_welch;
    OUT(iD).SnR_after_welch_all = SnR_after_welch_all;
    OUT(iD).SnR_r_val_after = SnR_r_val_after;
    OUT(iD).pwavelet_before = mean(Pwv,2);
    OUT(iD).pwavelet_before_no_ssr = mean(Pwv_no_ssr,2);
    OUT(iD).pwavelet_after = mean(Pwv_on,2);

    OUT(iD).pwavelet_fq = fq;
    OUT(iD).welch_fq = welch_fq;
    OUT(iD).pwelch_before = PW_before;
    OUT(iD).pwelch_after_all = PW_after;
    OUT(iD).pwelch_before_no_SSR = PW_before_no_SSR; % Use this as the template
    OUT(iD).Frequency_for_SigToNoise = Frequency_for_SigToNoise_BLK;
    OUT(iD).Peak_mV_from_hilbert_before = mean(abs(hilbert(LFP_off)));
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the frequency response during baseline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if PLOT_IT
        figure(iD*10)
        clf
        subplot(4,1,1)
        plot(LFP_off_x_sec, LFP_off)
        axis tight
        plot_horiz_line_at_zero(OUT(iD).Peak_mV_from_hilbert_before)
        
        xlabel('s'); ylabel('uV')
        pubify_figure_axis
        title([num2str(don(iD).Freq) ' Hz Baseline'])
        
        subplot(4,1,2)
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
        subplot(4,1,3)
        
        imagesc([],fq,Pwv); axis xy

        
        
        subplot(4,1,4)
        %plot(S_fq, mean(S_before,2))
        yyaxis left
        plot(welch_fq, 10*log10(PW_before))
        yyaxis right
        plot(fq,OUT(iD).pwavelet_before )
        hold on
        plot(fq,OUT(iD).pwavelet_after,'k')
        
        xlabel('Frequency (Hz)');ylabel('10*log10(squared uV/Hz)')
        axis tight
        title(sprintf('Welch PSD, Max Fq = %2.4f, SnR = %2.4f SnR Border Fq %2.1f %2.1f Hz ',pq_fq, SnR_before_welch, Frequency_for_SigToNoise_BLK(2), Frequency_for_SigToNoise_BLK(3)))
        pubify_figure_axis
        label_figure(mname);
        % Save figure.
        saveas(gcf,fullfile('Figures',sprintf('Before_Fq_%d.png',don(iD).Freq)),'png')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the ratio of the power in the outer frequencies compared to the
    % target frequency for each band.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First, iterate through each interval and determine the spectrogram
    % for each inter-pulse period.
    %     sc = [don(iD).EventRecs_sec(:) don(iD).EventRecsDown_sec(:)];
    down = don(iD).EventRecsDown_sec(:);
    %     dur_s = median(sc(2:end,1) - sc(1:end-1,1));
    %     newx = linspace(0, dur_s,round(dur_s*LFP_sFreq));
    samples_before = LFP_sFreq*.01;
    if FOR_HIGH_FREQUENCIES
        samples_after = LFP_sFreq*.180;
        skip=1;
    else
        % For slower frequencies.
        samples_after = LFP_sFreq*2.50;
        skip = 10;
    end
    OUT(iD).skip = skip;

    if subtract_template
        [L, ix, newx] = PETH_EEG_simple([don(iD).x_LFP_sec(:) don(iD).LFPsub(:)], ...
            down(1:skip:end), samples_before, samples_after,LFP_sFreq, false);
    else
        [L, ix, newx] = PETH_EEG_simple([don(iD).x_LFP_sec(:) don(iD).LFP(:)], ...
            down(1:skip:end), samples_before, samples_after,LFP_sFreq, false);
    end
    % get rid of the last trial
    L = L(5:end-5,:)';
    if Cols(L) > 10
        L = L(:,1:10); % Really - we don't need that many.
    end
    if any(isnan(L(:)))
        error('a')
    end
    welch_IX = newx > META.Time_from_start_to_ignore_sec; % & newx < newx(end)-0.01;
%     spec_IX  = newx > META.Time_from_start_to_ignore_sec; % & newx < newx(end)-0.01;
    
    %     function [M, ix, x] = PETH_EEG_simple(EEG_t_data, alignments_t, samples_before, samples_after,sFreq, PLOT_IT)
    %     L = zeros(length(newx), Rows(sc)-1);
    %     last_ix = round(length(newx)-(0.05*LFP_sFreq):length(newx));
    %     for iScan = 1:Rows(sc)-1
    %         IX = x_sec_after >= EventRecsEnds_sec_block(iScan) & x_sec_after <= EventRecs_sec_block(iScan+1);
    %         x = x_LFP_sec(IX) - EventRecs_sec_block(iScan);
    %         L(:,iScan) = interp1(x,LFP_after(IX),newx);
    %         L(~spec_IX,:) = nan;
    %         NANIX = ~isnan(L(:,iScan));
    %         L(NANIX,iScan) = detrend(L(NANIX,iScan));
    %         L(:,iScan) = L(:,iScan) - nanmean(L(last_ix,iScan));
    %     end
    if iD > 1
        b = fir1(95,Frequency_for_SigToNoise_BLK(2:3)*2/LFP_sFreq);
    else
        b = [];
    end
    %     [aa, bb ]= butter(5,Frequency_for_SigToNoise(iD,2:3)*2/LFP_sFreq);
    
    if 0
        b = fir1(148,Frequency_for_SigToNoise_BLK(2:3)*2/LFP_sFreq);
        [aa, bb ]= butter(5,Frequency_for_SigToNoise_BLK(2:3)*2/LFP_sFreq);
        figure
        freqz(b,1,512,LFP_sFreq)
        figure
        freqz(bb,aa,512,LFP_sFreq)
        LFP_low = filtfilt(b,1,LFP(ix_before));
    end
    % We need good time resolution...
    %     window_sec = min([5/block_fq 0.05]); %  a miniumum of around 5 cycles seems like a bare minium
   
    %     [~,fq,~,P] = spectrogram(LFP_off,window_samples,round(window_samples/3),SnR_FQ,LFP_sFreq);
    %     S_before = abs(P); %
    %     %     S_before = 10*log10(abs(P)); % From the matlab docs.
    %     x_S_before_sec = linspace(x_LFP_sec(ix_before(1)),x_LFP_sec(ix_before(end)),Cols(S_before));
    %     SnR_before = S_before(2,:)./mean(S_before([1 3],:));
    %     SNR_by_time = [];
    %     Peak_mV_from_hilbert_after = zeros(Cols(L),1);
    %     Peak_mV_from_hilbert_after2 = zeros(Cols(L),1);
    Lv = zeros(size(L));
    Lvf = zeros(size(L));
    
    window_samples_scan = ceil(LFP_sFreq * 3/OUT(iD).Target_Fq); % Make sure there are at least 3 cycles.
    
    for iScan = 1:Cols(L)
        %         v = detrend(L(welch_IX,iScan));
        v = L(welch_IX,iScan);
        if iD > 1
            v_filt1 = filtfilt(b,1,v);
        else
            v_filt1 = v;
        end
        
        %v_filt1 = filfilt(bb,aa,v);
        mid = round(length(v)/2);
        Lv(welch_IX,iScan) = v;
        Lvf(welch_IX,iScan) = v_filt1;
        
        if any(isnan(v_filt1))
            error('as')
        end
        
        [~ , ~, ~, P]  = spectrogram( v, window_samples,...
            round( window_samples/2), SnR_FQ, LFP_sFreq);
        
        if iScan == 1
            [~ , tfq, ~, FULP]  = spectrogram( v, window_samples,...
                round( window_samples/2), welch_fq, LFP_sFreq);
            fullP = zeros(size(FULP,1),size(FULP,2),Cols(L))*nan;
            fullP(:,:,iScan) = FULP;
            rmatch = NaN(Cols(L),Cols(FULP));
        else
            [~ , tfq, ~, FULP]  = spectrogram( v, window_samples,...
                round( window_samples/2), welch_fq, LFP_sFreq);
            fullP(:,:,iScan) = FULP;
        end
        for ii = 1:Cols(FULP)
            tmp = corrcoef(FULP(:,ii),OUT(iD).pwelch_before_no_SSR(:));
            rmatch(iScan,ii) = tmp(2);
        end
        if iScan == 3
            [~,Pwv] = SPEC_waveletdecomp(OUT(iD).pwavelet_fq,v,LFP_sFreq,wavelet_size_in_scan);
            [~,Pwvls] = SPEC_waveletdecomp(META.fq_logspace,v,LFP_sFreq,wavelet_size_in_scan);
            imagesc([],[], log10(Pwvls))
            set(gca,'YTick',1:length(META.fq_logspace),'YTickLabel',META.fq_logspace)
            
            %             Pwv = Pwv(:,round(2*LFP_sFreq):(end-round(2*LFP_sFreq)));
            for ii = 1:Cols(Pwv)
                tmp = corrcoef(Pwv(:,ii),OUT(iD).pwavelet_before_no_ssr(:));
                rmatch2(ii) = tmp(2);
            end
        end
        [PW] = pwelch( v, length(v), [], welch_fq, LFP_sFreq); % plot(welch_fq,PW)
        PWn = PW - PW_before;
        PWlog = 10*log10(PW);
        [PWsnr, ~] = pwelch(v, length(v), [],SnR_FQ , LFP_sFreq);
        PWsnr_n = PWsnr - PW_snr_before;
        
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
    
    if 1
        
        figure(20201)
        clf
        subplot(2,1,1)
        plot(mean(L'))
        hold on
        plot(mean(Lv'))
        plot(mean(Lvf'))
        legend('L','Lv','Lvf')
        subplot(2,1,2)
        plot_confidence_intervals(rmatch)
        title('rmatch')
    end
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
    
    OUT(iD).META = META;
    OUT(iD).LFP_sFreq = LFP_sFreq;
    OUT(iD).Pwvls = single(Pwvls);
    OUT(iD).rmatch = rmatch;
    OUT(iD).rmatch_mn = mean(rmatch);
    OUT(iD).rmatch_wv_mn = rmatch2;

    
    OUT(iD).Mean_Welch_Power_after = nanmean(pwelch_after_snr,2)';
    OUT(iD).Mean_Welch_Power_after_norm = nanmean(pwelch_after_snr_norm,2)';
    OUT(iD).Mean_Welch_SnR_after = nanmean(SnR_after_welch);
    OUT(iD).Mean_Welch_SnR_after_norm = nanmean(SnR_after_welch_norm);
    OUT(iD).Peak_Fq_after = pq_fq;
    OUT(iD).Mean_Full_Welch_after = mean(pwelch_after,2);
    OUT(iD).SNR_by_time =mean(SNR_by_time,2);
    OUT(iD).SNR_by_time_x_sec =newx(range_to_test);
    OUT(iD).Peak_mV_from_hilbert_after = median(Peak_mV_from_hilbert_after);
    %%%%%%%%%%%%%%%%%%%%
    if 0
        
        figure(iD*100)
        clf
        subplot(3,2,1)
        plot(newx,L)
        axis tight; xlabel('s');ylabel('uV'); pubify_figure_axis
        hold on
        plot(0,OUT(iD).Peak_mV_from_hilbert_after,'r>');
%         a = axis; a(1:2) = [0 dur_s]; axis(a);
        title([num2str(block_fq) ' Hz'])
        
        pubify_figure_axis
        subplot(3,2,2)
        
        plot(x_S_after_sec,S_after,'o-')
        axis tight; xlabel('s');ylabel('Power'); pubify_figure_axis
        title(sprintf('Max Fq = %f Hz',pq_fq))
        legend('lower','target','upper')
        
        subplot(3,2,3)
        % Spectrogram.
        cla
        x = linspace(x_S_after_sec(1) + META.Time_from_start_to_ignore_sec,...
            x_S_after_sec(end),size(fullP,2));
        V = squeeze(nanmean(real(fullP),3));
        imagesc(x, tfq,V)
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
        saveas(gcf,fullfile(data_dir,'Figures',sprintf('Fq %f_summary.png',block_fq)),'png')
        
    end
    fprintf('.')
end
% Save pre data.
disp('Saving')
save(mname,'OUT');

%%
R = [];
for ii = 1:length(OUT)
      R(ii,:) =  OUT(ii).rmatch_wv_mn;
%      R(ii,:) =  OUT(ii).rmatch_mn;
end

SnR_no_FSCV = [OUT.SnR_before_welch];
SnR = [OUT.Mean_Welch_SnR_after];
FQ  = [OUT.Target_Fq];
xIX = FQ < 100;
ix = 1:LFP_sFreq/5;

figure
imagesc((1:length(ix))/LFP_sFreq,[],R(:,ix).^2)
set(gca,'YTick',1:length(FQ))
set(gca,'YTickLabel',FQ)
colorbar_label('r^2')
colormap(jet);
figure
plot((1:length(ix))/LFP_sFreq,R(:,ix)'.^2)

%%
clr = lines(2);
figure
bar(log(SnR),'FaceColor','k')
hold on
tmp = SnR_no_FSCV;
tmp(end+1) = tmp(end);
stairs((1:length(tmp))-.5,log(tmp),'Color',clr(2,:))
set(gca,'XTick',1:length(SnR))
set(gca,'XTickLabel',FQ)
set(gca,'XTickLabelRotation',45)
a = axis;
a(3) = -1;
axis(a)
ylabel('log(SnR)')
xlabel('Hz')
legend('Snr','Snr no FSCV')
grid on
plot_horiz_line_at_zero
label_figure(mname)

%%
figure
bar(SnR(xIX),'k')
set(gca,'XTick',1:length(SnR(xIX)))
set(gca,'XTickLabel',FQ(xIX))
plot_horiz_line_at_zero();
pubify_figure_axis
ylabel('SnR')
xlabel('Frequency (Hz)')
grid on
label_figure(mname);

