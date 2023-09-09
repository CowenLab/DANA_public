function [OUT,META] = DANA_analyze_brain_gel_post_4_14_16_v3()
% Create brain gel figures for the DANA project.
%  - assumes INTAN support functions are in the path.
% Assumes that you are in a directory with the rhd, dat, and datlfp data.
% or pass in the data_dir with the .dat files:
% e.g., data_dir = 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF on';
%
%
% Cowen 2016 (with some inspiration from JP's code)d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'defaultAxesFontName', 'Arial');
mname = 'DANA_analyze_brain_gel_post_4_14_16_v3';

% savetype = 'png';
savetype = 'fig';
saveit = false;
PLOT_IT = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( 'C:\Cowen\Data\DANA\04-14-2016 - brain gel\no SSR WF on','dir');
    root  = 'C:\Cowen\Data\DANA\';
else
    root  = 'G:\Cowen\Data\DANA\';
end
data_dir_wfon_no_ssr = [ root '04-14-2016 - brain gel\no SSR WF on'];
data_dir_wfoff_no_ssr = [ root '04-14-2016 - brain gel\no SSR WF off'];
data_dir_wfon_ssr = [ root '04-14-2016 - brain gel\with SSR WF on'];
data_dir_wfoff_ssr = [ root '04-14-2016 - brain gel\with SSR  WF off'];
% do not do this right now. Not so great with slow frequency artifact. Nice
% for HF noise though. Maybe post-filtering.
%%
wavelet_size = 5; % 10 or 20 works as well. Higher and get better fq resolution. Log so maybe not as importnt.
wavelet_size_in_scan = wavelet_size; % These two need to be the same for now. Problem - values probably need to change radically depending on the frequency - lower for low frequencies, higher for high frequencies.

% pad_future_sec = 0.005;
pad_future_sec = 0.005;
% META.Time_from_start_to_ignore_sec = 0.01;
META.Time_from_start_to_ignore_sec = pad_future_sec;
META.fq_logspace = unique((logspace(log10(1),log10(3000),1650)));
META.fq_logspace(META.fq_logspace < 1) = [];
META.fq_low_range = 1:.2:55;
META.fq_high_range_log = unique(ceil(logspace(log10(50),log10(3000),300)));
decimation_factor = 3; % set to 1 for NO decimation
OUT.META = META;
META.fq_HI = 600:20:6000;
META.fq_LO = 0.5:.5:250;

% META.post_scan_dur_sec = .191;
META.post_scan_dur_sec = .175; % Why this? Why not .191? To deal with slop perhaps.

data_dir = root;
figdir = fullfile(data_dir,'Figures');
if ~exist(fullfile(data_dir,'Figures'),'dir')
    mkdir(fullfile(data_dir,'Figures'))
end
%% Create a list of sub directories.
fnames = {'amp-D-000.dat' 'amp-D-001.dat' 'amp-D-002.dat' 'amp-D-003.dat' 'amp-D-012.dat' 'amp-D-013.dat' 'amp-D-014.dat' 'amp-D-015.dat' };
% fnames = {'amp-D-000.dat'};
for iFile = 1:length(fnames)
    don_ssr = DANA_Load_brain_gel_data_dirs(data_dir_wfon_ssr, fnames{iFile}, 'board-DIN-03.dat',decimation_factor,false);
    doff_ssr = DANA_Load_brain_gel_data_dirs(data_dir_wfoff_ssr, fnames{iFile}, 'board-DIN-03.dat',decimation_factor,false);
    don_nossr = DANA_Load_brain_gel_data_dirs(data_dir_wfon_no_ssr, fnames{iFile}, 'board-DIN-03.dat',decimation_factor,false);
    doff_nossr = DANA_Load_brain_gel_data_dirs(data_dir_wfoff_no_ssr, fnames{iFile}, 'board-DIN-03.dat',decimation_factor,false);
    %%
    FQ = [don_nossr.Freq];
    % make sure that the psd frequencies include the specific targets...
    META.fq_logspace = unique([META.fq_logspace FQ]);
    META.fq_high_range_log = unique([META.fq_high_range_log FQ]);
    META.fq_high_range_log(META.fq_high_range_log < 50) = [];
    %%
    
    %           for iFreq = 6 %1:length(don_ssr) % 4 = 5hz ,5 = 7 hz.
    for iFreq = 1:length(don_ssr)
        Target_Fq = don_nossr(iFreq).Freq;
        OUT(iFreq).Target_Fq = Target_Fq;
        OUT(iFreq).LFP_sFreq = don_nossr(iFreq).LFP_sFreq;
        
        tit = [ num2str(Target_Fq) ' Hz'];
        % Align recs to the stim - this will make plots nicer.
        strec = 10;
        don_ssr(iFreq).x_LFP_sec = don_ssr(iFreq).x_LFP_sec - don_ssr(iFreq).EventRecs_sec(strec);
        don_ssr(iFreq).EventRecsDown_sec = don_ssr(iFreq).EventRecsDown_sec - don_ssr(iFreq).EventRecs_sec(strec);
        don_ssr(iFreq).EventRecs_sec = don_ssr(iFreq).EventRecs_sec - don_ssr(iFreq).EventRecs_sec(strec);
        doff_ssr(iFreq).x_LFP_sec = doff_ssr(iFreq).x_LFP_sec - don_ssr(iFreq).EventRecs_sec(strec);
        don_nossr(iFreq).x_LFP_sec = don_nossr(iFreq).x_LFP_sec - don_nossr(iFreq).EventRecs_sec(strec);
        don_nossr(iFreq).EventRecsDown_sec = don_nossr(iFreq).EventRecsDown_sec - don_nossr(iFreq).EventRecs_sec(strec);
        don_nossr(iFreq).EventRecs_sec = don_nossr(iFreq).EventRecs_sec - don_nossr(iFreq).EventRecs_sec(strec);
        doff_nossr(iFreq).x_LFP_sec = doff_nossr(iFreq).x_LFP_sec - don_nossr(iFreq).EventRecs_sec(strec);
        
        % Convert everything to mV.
        don_ssr(iFreq).LFP = don_ssr(iFreq).LFP/1000;
        doff_ssr(iFreq).LFP = doff_ssr(iFreq).LFP/1000;
        don_nossr(iFreq).LFP = don_nossr(iFreq).LFP/1000;
        doff_nossr(iFreq).LFP = doff_nossr(iFreq).LFP/1000;
        
        LFP_sFreq = don_ssr(iFreq).LFP_sFreq;
        step_size_samples = LFP_sFreq/200; % so 5 ms bin size for steps
        
        bpFilt = designfilt('bandpassiir', 'FilterOrder', 12, 'PassbandFrequency1', 100, 'PassbandFrequency2', 160, 'PassbandRipple', 1, 'SampleRate', LFP_sFreq);
        don_nossr(iFreq).LFP_filt140 = filtfilt(bpFilt,double(don_nossr(iFreq).LFP ));
        don_ssr(iFreq).LFP_filt140 = filtfilt(bpFilt,double(don_ssr(iFreq).LFP ));
        
        IX_ON = don_nossr(iFreq).x_LFP_sec >= 5 & don_nossr(iFreq).x_LFP_sec <= 5.5;
        IX_OFF = doff_nossr(iFreq).x_LFP_sec >= 5 & doff_nossr(iFreq).x_LFP_sec <= 5.5;
        
        figure(1)
        clf
        subplot(2,1,1)
        plot(don_nossr(iFreq).x_LFP_sec(IX_ON), don_nossr(iFreq).LFP(IX_ON),'LineWidth',3)
        hold on
        plot(doff_nossr(iFreq).x_LFP_sec(IX_OFF), doff_nossr(iFreq).LFP(IX_OFF),'LineWidth',3)
        legend('WF On','WF Off'); legend boxoff
        title([ 'No SSR ' tit ])
        axis tight
        
        pubify_figure_axis
        %         ylabel('$\mu$V','Interpreter','Latex')
        ylabel('mV','Interpreter','Latex')
        xlabel('sec')
        
        
        IX_ON = don_ssr(iFreq).x_LFP_sec >= 5 & don_ssr(iFreq).x_LFP_sec <= 5.5;
        IX_ONlong = don_ssr(iFreq).x_LFP_sec >= 5 & don_ssr(iFreq).x_LFP_sec <= 10;
        
        IX_OFF = doff_ssr(iFreq).x_LFP_sec >= 5 & doff_ssr(iFreq).x_LFP_sec <= 5.5;
        
        subplot(2,1,2)
        plot(don_ssr(iFreq).x_LFP_sec(IX_ON), don_ssr(iFreq).LFP(IX_ON),'LineWidth',3)
        hold on
        plot(doff_ssr(iFreq).x_LFP_sec(IX_OFF), doff_ssr(iFreq).LFP(IX_OFF),'LineWidth',3)
        legend('WF On','WF Off'); legend boxoff
        axis tight
        title([ 'SSR ' tit ])
        pubify_figure_axis
        ylabel('mV','Interpreter','Latex')
        xlabel('sec')
        
        %         label_figure(mname)
        if saveit, saveas(gcf,fullfile(figdir,sprintf('SinglesubplotExample%s',tit)),savetype), end
        %%%
        
        figure(12)
        clf
        [OUT(iFreq).welchssr, OUT(iFreq).welchnossr, OUT(iFreq).welchpure] = DANA_plot_psd_ssr_on_off(don_ssr(iFreq).LFP(IX_ONlong), don_nossr(iFreq).LFP(IX_ONlong), META.fq_logspace, LFP_sFreq,Target_Fq,doff_nossr(iFreq).LFP(IX_ONlong));
        if saveit, saveas(gcf,fullfile(figdir,sprintf('welchssrnossr%s',tit)),savetype), end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % plot template match
        % but choose place so that it starts with offset.
        
        ix = find(don_nossr(iFreq).EventRecsDown_sec > 5,1,'first');
        st_sec = don_nossr(iFreq).EventRecsDown_sec(ix);
        IX_ON = don_nossr(iFreq).x_LFP_sec >= st_sec & don_nossr(iFreq).x_LFP_sec <= st_sec + META.post_scan_dur_sec;
        IX_OFF = doff_nossr(iFreq).x_LFP_sec >= st_sec & doff_nossr(iFreq).x_LFP_sec <= st_sec + META.post_scan_dur_sec;
        
        [LFP_target] = Artificial_LFP(LFP_sFreq, 2, Target_Fq, 0, 0 );
        LFP_target = LFP_target*rms(don_nossr(iFreq).LFP(IX_OFF));
        
        %         figure(1)
        %         [ OUT(iFreq).rmatch_no_ssr ,OUT(iFreq).rmatch_ideal,INFO]= SPEC_time_to_recover_signal(don_nossr(iFreq).LFP(IX_ON), LFP_sFreq, LFP_target,META.fq_logspace,step_size_samples,'wavelet',true);
        %         OUT(iFreq).rmatch_no_ssr_x_sec = INFO.steps/LFP_sFreq;
        %
        OUT(iFreq).rmatch_nossr_HF = [];
        OUT(iFreq).rmatch_nossr_LF = [];
        
        if Target_Fq > 400
            [OUT(iFreq).rmatch_nossr_HF,OUT(iFreq).rmatch_nossr_ideal_HF,INFO]= SPEC_time_to_recover_signal(don_nossr(iFreq).LFP(IX_ON), LFP_sFreq, LFP_target,META.fq_HI,step_size_samples,'wavelet',true);
            OUT(iFreq).rmatch_no_ssr_ttr_sec_HF = INFO.time_to_recover_s;
            OUT(iFreq).rmatch_no_ssr_ttr_ideal_sec_HF = INFO.time_to_recover_ideal_s;
        end
        if Target_Fq < 300
            [OUT(iFreq).rmatch_nossr_LF,OUT(iFreq).rmatch_nossr_ideal_LF,INFO]= SPEC_time_to_recover_signal(don_nossr(iFreq).LFP(IX_ON), LFP_sFreq, LFP_target,META.fq_LO,step_size_samples,'wavelet',true);
            OUT(iFreq).rmatch_no_ssr_ttr_sec_LF = INFO.time_to_recover_s;
            OUT(iFreq).rmatch_no_ssr_ttr_ideal_sec_LF = INFO.time_to_recover_ideal_s;
        end
        if Target_Fq >20
            disp('stpo')
        end
        %
        OUT(iFreq).rmatch_no_ssr_x_sec = INFO.x_st;
        %
        %         figure(2)
        %         [OUT(iFreq).rmatch_pure,~,~,~,~,steps]= SPEC_time_to_recover_signal(doff_nossr(iFreq).LFP(IX_OFF), LFP_sFreq, LFP_target,META.fq_logspace,step_size_samples,true);
        %         OUT(iFreq).rmatch_x_sec_pure = steps/LFP_sFreq;
        OUT(iFreq).rmatch_pure_HF = [];
        if Target_Fq > 400
            [OUT(iFreq).rmatch_pure_HF,~,INFO]= SPEC_time_to_recover_signal(don_nossr(iFreq).LFP(IX_OFF), LFP_sFreq, LFP_target,META.fq_HI,step_size_samples,'wavelet',false);
            OUT(iFreq).rmatch_pure_ttr_sec_HF = INFO.time_to_recover_s;
            OUT(iFreq).rmatch_pure_ttr_ideal_sec_HF = INFO.time_to_recover_ideal_s;
        end
        OUT(iFreq).rmatch_pure_LF = [];
        if Target_Fq < 300
            [OUT(iFreq).rmatch_pure_LF,~,INFO]= SPEC_time_to_recover_signal(don_nossr(iFreq).LFP(IX_OFF), LFP_sFreq, LFP_target,META.fq_LO,step_size_samples,'wavelet',false);
            OUT(iFreq).rmatch_pure_ttr_sec_LF = INFO.time_to_recover_s;
            OUT(iFreq).rmatch_pure_ttr_ideal_sec_LF = INFO.time_to_recover_ideal_s;
        end
        OUT(iFreq).rmatch_no_ssr_x_sec = INFO.x_st;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % plot template match
        % convert stim periods to nan.
        LFP = don_ssr(iFreq).LFP;
        st = don_ssr(iFreq).EventRecs_sec - 0.001;
        ed = don_ssr(iFreq).EventRecsDown_sec + pad_future_sec;
        st = st(1:length(ed));
        for ii = 1:length(st)
            IX = don_ssr(iFreq).x_LFP_sec >= st(ii) &  don_ssr(iFreq).x_LFP_sec < ed(ii);
            LFP(IX) = 0;
        end
         figure
         plot( don_ssr(iFreq).LFP,'r');
         hold on
         plot(LFP,'b')
        
        ix = find(don_ssr(iFreq).EventRecsDown_sec > 5,1,'first');
        st_sec = don_ssr(iFreq).EventRecsDown_sec(ix) + pad_future_sec;
        IX_ON = don_ssr(iFreq).x_LFP_sec >= st_sec & don_ssr(iFreq).x_LFP_sec <= st_sec + META.post_scan_dur_sec;
        %         IX_OFF = doff_ssr(iFreq).x_LFP_sec >= st_sec & doff_ssr(iFreq).x_LFP_sec <= st_sec + 1;
        
        [LFP_target] = Artificial_LFP(LFP_sFreq, 2, Target_Fq, 0, 0 );
        LFP_target = LFP_target*.5;
        
        %         figure(4)
        %         %             [rmatch_ssr,rmatch_ssr_ideal,M,Mideal,template,steps]= SPEC_time_to_recover_signal(don_ssr(iFreq).LFP(IX_ON), LFP_sFreq, LFP_target,META.fq_logspace,step_size_samples,true);
        %         [OUT(iFreq).rmatch_ssr,~,~,~,~,steps]= SPEC_time_to_recover_signal(LFP(IX_ON), LFP_sFreq, LFP_target,META.fq_logspace,step_size_samples,true);
        %         OUT(iFreq).rmatch_x_sec_ssr = steps/LFP_sFreq;
        
        OUT(iFreq).rmatch_ssr_HF  = [];
        if Target_Fq > 400
            [OUT(iFreq).rmatch_ssr_HF,~,INFO]= SPEC_time_to_recover_signal(detrend(LFP(IX_ON)), LFP_sFreq, ...
                LFP_target,META.fq_HI,step_size_samples,'wavelet',false);
            OUT(iFreq).rmatch_ssr_ttr_sec_HF = INFO.time_to_recover_s;
            OUT(iFreq).rmatch_ssr_ttr_ideal_sec_HF = INFO.time_to_recover_ideal_s;
        end
        OUT(iFreq).rmatch_ssr_LF = [];
        if Target_Fq < 300
            ix = find(don_ssr(iFreq).EventRecsDown_sec > 5);
            ttr = [];ttri = [];
            for jj = 1:10 % may be useless - artifact seemed pretty regular.
                st_sec = don_ssr(iFreq).EventRecsDown_sec(ix(jj)) + pad_future_sec;
                IX_ON = don_ssr(iFreq).x_LFP_sec >= st_sec & don_ssr(iFreq).x_LFP_sec <= st_sec + META.post_scan_dur_sec;

                [RR(jj,:),~,INFO]= SPEC_time_to_recover_signal(detrend(LFP(IX_ON)), LFP_sFreq, ...
                    LFP_target,META.fq_LO,step_size_samples,'wavelet',false);
                ttr(jj) = INFO.time_to_recover_s;
                ttri(jj) = INFO.time_to_recover_ideal_s;
            end
            OUT(iFreq).rmatch_ssr_LF = nanmean(RR);
            OUT(iFreq).rmatch_ssr_ttr_sec_LF = nanmean(ttr);
            OUT(iFreq).rmatch_ssr_ttr_ideal_sec_LF = nanmean(ttri);
        end
        OUT(iFreq).rmatch_ssr_x_sec = INFO.x_st;

        % HF only
        
        OUT(iFreq).Target_Fq = Target_Fq;
        OUT(iFreq).LFP_sFreq = LFP_sFreq;
        drawnow
        
        fprintf('.')
    end
    % Save pre data.
    disp('Saving')
    [~,fn] = fileparts(fnames{iFile});
    save([mname fn],'OUT','mname');
end
% load(mname)
%
%
%% Template matching SnR analysis.
savetype = 'eps';
FQ  = [don_ssr.Freq];
fIX = FQ >= 1;
LFP_sFreq = don_ssr(1).LFP_sFreq;
nCol = length(OUT(1).rmatch_ssr_LF);
RpureLF = zeros(length(OUT),nCol);
RssrHF = zeros(length(OUT),nCol);
Rnossr = zeros(length(OUT),nCol);
RidealLF = zeros(length(OUT),nCol);
RnossrHF = zeros(length(OUT),nCol);
RnossrLF = zeros(length(OUT),nCol);
RssrLF = zeros(length(OUT),nCol);
RpureHF = zeros(length(OUT),nCol);
for ii = 1:length(OUT)
    %     Rssr(ii,:) =  OUT(ii).rmatch_ssr';
    %     Rnossr(ii,:) =  OUT(ii).rmatch_no_ssr';
    %     Rpure(ii,:) =  OUT(ii).rmatch_pure';
    if ~isempty(OUT(ii).rmatch_ssr_HF)
        RssrHF(ii,:) =  OUT(ii).rmatch_ssr_HF';
    end
    if ~isempty(OUT(ii).rmatch_ssr_LF)
        RssrLF(ii,:) =  OUT(ii).rmatch_ssr_LF';
    end
    if ~isempty(OUT(ii).rmatch_nossr_HF)
        RnossrHF(ii,:) =  OUT(ii).rmatch_nossr_HF';
    end
    if ~isempty(OUT(ii).rmatch_nossr_LF)
        RnossrLF(ii,:) =  OUT(ii).rmatch_nossr_LF';
    end
%     if ~isempty(OUT(ii).rmatch_pure_HF)
%         RpureHF(ii,:) =  OUT(ii).rmatch_pure_HF';
%     end
    if ~isempty(OUT(ii).rmatch_pure_LF)
        RpureLF(ii,:) =  OUT(ii).rmatch_pure_LF';
    end
    if ~isempty(OUT(ii).rmatch_nossr_ideal_LF)
        RidealLF(ii,:) =  OUT(ii).rmatch_nossr_ideal_LF';
    end
    
end
RssrHF = RssrHF(fIX,:);
RnossrHF = RnossrHF(fIX,:);
RnossrLF = RnossrLF(fIX,:);
% RpureHF = RpureHF(fIX,:);
RssrLF = RssrLF(fIX,:);
RidealLF = RidealLF(fIX,:);
RpureLF = RpureLF(fIX,:);

%%
x_s = OUT(ii).rmatch_no_ssr_x_sec;
newFQ = FQ(fIX);
% ix = 1:LFP_sFreq/5;
ix = 1:length(x_s);
flab = [];
for ii = 1:length(newFQ)
    flab{ii} = [num2str(newFQ(ii)) 'Hz'];
end

subplot(2,3,1)
imagesc(x_s(ix),[],RnossrHF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('NOSSR HF')
% set(gca,'XLim',[0 .18])


subplot(2,3,2)
imagesc(x_s(ix),[],RnossrLF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('NO SSR LF')
% set(gca,'XLim',[0 .18])


subplot(2,3,3)
imagesc(x_s(ix),[],RidealLF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('RidealLF LF')
% set(gca,'XLim',[0 .18])



subplot(2,3,4)
imagesc(x_s(ix),[],RssrHF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('SSR HF')
% set(gca,'XLim',[0 .18])


subplot(2,3,5)
imagesc(x_s(ix),[],RssrLF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('SSR LF')
% set(gca,'XLim',[0 .18])

subplot(2,3,6)
imagesc(x_s(ix),[],RpureLF(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('RpureLF')
% set(gca,'XLim',[0 .18])



ftyp = {'png' 'eps' 'fig'};
for ii = 1:length(ftyp)
    saveas(gcf,sprintf('./Figures/FtemplateSnRlines.%s',ftyp{ii}))
end
%%
figure
% subplot(1,2,1)

set(0,'DefaultAxesColorOrder',magma(length(newFQ)))

plot((1:length(ix))/LFP_sFreq,Rnossr(:,ix)'.^2,'Linewidth',2)
% label_figure(mname);
% label_figure(fname,'bottom left');
legend(flab)
legend boxoff
xlabel('Seconds')
ylabel('r^2')
txt = '';
title(txt)

saveas(gcf,fullfile(pwd,['FtemplateSnRlines ' txt]),savetype)


%% Give Kate the vector of uV and time for the 0Hz SSR...


% y_uVsq = OUT(1).welchnossr;
% y_uVsq = OUT(1).welchpure;
figure

subplot(1,2,1)
x_Hz = OUT(1).wavelet_psd_xlabs;
y_uVsq = OUT(1).welchssr*1000;
plot(x_Hz,y_uVsq)
set(gca,'XLim',[0 150])
xlabel('Hz')
ylabel('uV^2')
title('SSR')
pubify_figure_axis
saveas(gcf,'ssr_no_input.eps')
csvwrite('ssr_out.csv',[x_Hz(:) y_uVsq(:)]);
subplot(1,2,2)
x_Hz = OUT(1).wavelet_psd_xlabs;
y_uVsq = OUT(1).welchpure*1000;;
plot(x_Hz,y_uVsq)
set(gca,'XLim',[0 150])
xlabel('Hz')
ylabel('uV^2')
title('NO SSR')
pubify_figure_axis
saveas(gcf,'nossr_no_input.eps')
csvwrite('nossr_out.csv',[x_Hz(:) y_uVsq(:)]);







%%

fqs_to_use = [1, 5 , 7, 50, 100, 1000];
FQ = [OUT.Target_Fq];
[f,fqix] = intersect(FQ,fqs_to_use);
nc = ceil(length(fqs_to_use)/2);
h = [];
logscale = false;

figure(1010)
clf
for ii = 1:length(fqs_to_use)
    fq = fqs_to_use(ii);
    fq_ix = fqix(ii);
    x = OUT(fq_ix).wavelet_psd_xlabs;
    h(ii) = subplot(2,nc,ii);
    plot( x, (OUT(fq_ix).welchpure)','g','LineWidth',4);
    hold on
    plot( x, (OUT(fq_ix).welchssr)','b','LineWidth',2);
    plot( x, (OUT(fq_ix).welchnossr)','k','LineWidth',2);
    
    if OUT(fq_ix).Target_Fq <  50
        xl = [0 20];
    elseif  OUT(fq_ix).Target_Fq >=  40 && OUT(fq_ix).Target_Fq < 80
        xl = [20 80];
    elseif  OUT(fq_ix).Target_Fq >  80 && OUT(fq_ix).Target_Fq < 200
        xl = [60 150];
    elseif  OUT(fq_ix).Target_Fq >=  200
        xl = [900 1100];
    end
    axis tight
    
    if logscale
        set(gca,'YScale','log')
        ylabel('$$ log10 \mu V^2 $$','Interpreter','latex')
        txt = 'log';
    else
        ylabel('$$ \mu V^2 $$','Interpreter','latex')
        txt = '';
    end
    
    a = axis;
    
    %     plot(fq,a(4)*1.1,'rv','MarkerSize',8,'MarkerFaceColor','r')
    plot(fq,a(3),'r^','MarkerSize',8,'MarkerFaceColor','r')
    
    
    set(gca,'XLim',xl)
    if ii == 1
        legend('pure','ssr','nossr')
        %         legend('pure','nossr')
        legend boxoff
    end
    title([num2str(OUT(fq_ix).Target_Fq) ' Hz'])
    xlabel('Hz')
end
% max([get(h,'Ylim')])
equalize_y_axes(h)
% label_figure(mname);
ftyp = {'png' 'eps' 'fig'};
for ii = 1:length(ftyp)
    saveas(gcf,sprintf('./Figures/Welch_spectraNOSSR%s.%s',txt,ftyp{ii}))
end

%% Now we need the time to recover plot.
R = [];
fq_ix = 1:14;
HF_ONLY = false;
% fq_ix = 1:17;
flabs = flab(fq_ix);
R{1} = Rssr(fq_ix,:).^2;
R{2} = Rnossr(fq_ix,:).^2;
R{3} = Rpure(fq_ix,:).^2;
%   x_s = OUT(1).rmatch_x_sec_ssr(1:end-1);
x_s = OUT(1).rmatch_ssr_high_fq_x_sec(1:end-1);
if HF_ONLY
    % Special case: Just look at High-frequencies and time to recover.
    fq_ix = 15:length(flab);
    flabs = flab(fq_ix);
    R = [];
    R{1} = RssrHF(fq_ix,:).^2;
    R{2} = RnossrHF(fq_ix,:).^2;
    R{3} = RpureHF(fq_ix,:).^2;
end


x_s = x_s - x_s(1);
th = 0.7;
clrs = hsv(length(R));
axis_meth = 2;
all_first = [];
figure(10101)
clf
for iR = 1:length(R)
    
    firsts = [];
    for ii = 1:Rows(R{iR})
        ix = find(R{iR}(ii,:) > th,1,'first');
        firsts(ii,1) = ii;
        if ~isempty(ix)
            %             plot(x_s(ix),ii,'r+')
            firsts(ii,2) = x_s(ix(1));
        else
            firsts(ii,2) = nan;
        end
    end
    %     firsts = firsts(~isnan(firsts(:,2)),:);
    switch axis_meth
        case 1
            plot(firsts(:,2), firsts(:,1),'.-','LineWidth',4,'MarkerSize',25,'Color',clrs(iR,:))
        case 2
            plot(firsts(:,1), firsts(:,2),'.-','LineWidth',4,'MarkerSize',25,'Color',clrs(iR,:))
    end
    hold on
    all_first(:,iR) = firsts(:,2);
    
end
%
csvwrite('time_to_recover_HF.csv',all_first)
% axis tight
pubify_figure_axis
legend('ssr','no ssr','pure')
legend boxoff
switch axis_meth
    case 1
        set(gca,'XLim',[0 .18])
        set(gca,'YTick',1:length(flabs))
        set(gca,'YTickLabel',flabs)
        xlabel('Seconds to Detection')
        ylabel('Target Frequency')
    case 2
        set(gca,'YLim',[0 .18])
        set(gca,'XTick',1:length(flabs))
        set(gca,'XTickLabel',flabs)
        ylabel('Seconds to Detection')
        xlabel('Target Frequency')
        set(gca,'XTickLabelRotation',45)
end

ftyp = {'png' 'eps' 'fig'};
for ii = 1:length(ftyp)
    saveas(gcf,sprintf('./Figures/Seconds_to_detection%d.%s',axis_meth,ftyp{ii}))
end