function [OUT,META] = DANA_analyze_brain_gel_post_4_14_16_v2()
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
mname = 'DANA_analyze_brain_gel_post_4_14_16_v2a';

% savetype = 'png';
savetype = 'fig';
saveit = false;
PLOT_IT = true;
PLOT_WAVELET = false;
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

META.Time_from_start_to_ignore_sec = 0.01;
META.Time_from_start_to_ignore_sec = 0.002;
META.fq_logspace = unique((logspace(log10(1),log10(3000),1650)));
META.fq_logspace(META.fq_logspace < 1) = [];
META.fq_low_range = 1:.2:55;
META.fq_high_range_log = unique(ceil(logspace(log10(50),log10(3000),300)));
decimation_factor = 3; % set to 1 for NO decimation
OUT.META = META;

data_dir = root;
figdir = fullfile(data_dir,'Figures');
if ~exist(fullfile(data_dir,'Figures'),'dir')
    mkdir(fullfile(data_dir,'Figures'))
end
%% Create a list of sub directories.
fnames = {'amp-D-000.dat' 'amp-D-001.dat' 'amp-D-002.dat' 'amp-D-003.dat' 'amp-D-012.dat' 'amp-D-013.dat' 'amp-D-014.dat' 'amp-D-015.dat' };
% fnames = {'amp-D-000.dat'};
for iF = 1:length(fnames)
    don_ssr = DANA_Load_brain_gel_data_dirs(data_dir_wfon_ssr, fnames{iF}, 'board-DIN-03.dat',decimation_factor,false);
    doff_ssr = DANA_Load_brain_gel_data_dirs(data_dir_wfoff_ssr, fnames{iF}, 'board-DIN-03.dat',decimation_factor,false);
    don_nossr = DANA_Load_brain_gel_data_dirs(data_dir_wfon_no_ssr, fnames{iF}, 'board-DIN-03.dat',decimation_factor,false);
    doff_nossr = DANA_Load_brain_gel_data_dirs(data_dir_wfoff_no_ssr, fnames{iF}, 'board-DIN-03.dat',decimation_factor,false);
    %%
    FQ = [don_nossr.Freq];
    % make sure that the psd frequencies include the specific targets...
    META.fq_logspace = unique([META.fq_logspace FQ]);
    META.fq_high_range_log = unique([META.fq_high_range_log FQ]);
    META.fq_high_range_log(META.fq_high_range_log < 50) = [];
    %%
    
%     for iD = 13 %1:length(don_ssr) % 4 = 5hz ,5 = 7 hz.
    for iD = 1:length(don_ssr)
        Target_Fq = don_nossr(iD).Freq;
        OUT(iD).Target_Fq = Target_Fq;
        OUT(iD).LFP_sFreq = don_nossr(iD).LFP_sFreq;
        
        tit = [ num2str(Target_Fq) ' Hz'];
        % Align recs to the stim - this will make plots nicer.
        strec = 10;
        don_ssr(iD).x_LFP_sec = don_ssr(iD).x_LFP_sec - don_ssr(iD).EventRecs_sec(strec);
        don_ssr(iD).EventRecsDown_sec = don_ssr(iD).EventRecsDown_sec - don_ssr(iD).EventRecs_sec(strec);
        don_ssr(iD).EventRecs_sec = don_ssr(iD).EventRecs_sec - don_ssr(iD).EventRecs_sec(strec);
        doff_ssr(iD).x_LFP_sec = doff_ssr(iD).x_LFP_sec - don_ssr(iD).EventRecs_sec(strec);
        don_nossr(iD).x_LFP_sec = don_nossr(iD).x_LFP_sec - don_nossr(iD).EventRecs_sec(strec);
        don_nossr(iD).EventRecsDown_sec = don_nossr(iD).EventRecsDown_sec - don_nossr(iD).EventRecs_sec(strec);
        don_nossr(iD).EventRecs_sec = don_nossr(iD).EventRecs_sec - don_nossr(iD).EventRecs_sec(strec);
        doff_nossr(iD).x_LFP_sec = doff_nossr(iD).x_LFP_sec - don_nossr(iD).EventRecs_sec(strec);
        
        % Convert everything to mV.
        don_ssr(iD).LFP = don_ssr(iD).LFP/1000;
        doff_ssr(iD).LFP = doff_ssr(iD).LFP/1000;
        don_nossr(iD).LFP = don_nossr(iD).LFP/1000;
        doff_nossr(iD).LFP = doff_nossr(iD).LFP/1000;
        
        LFP_sFreq = don_ssr(iD).LFP_sFreq;
        
        bpFilt = designfilt('bandpassiir', 'FilterOrder', 12, 'PassbandFrequency1', 100, 'PassbandFrequency2', 160, 'PassbandRipple', 1, 'SampleRate', LFP_sFreq);
        don_nossr(iD).LFP_filt140 = filtfilt(bpFilt,double(don_nossr(iD).LFP ));
        don_ssr(iD).LFP_filt140 = filtfilt(bpFilt,double(don_ssr(iD).LFP ));
        
        if 0
            
            plot(don_nossr(iD).LFP )
            hold on
            plot(don_nossr(iD).LFP_filt140 )
            plot(don_ssr(iD).LFP )
            plot(don_ssr(iD).LFP_filt140 )
        end
            IX_ON = don_nossr(iD).x_LFP_sec >= 5 & don_nossr(iD).x_LFP_sec <= 5.5;
            IX_OFF = doff_nossr(iD).x_LFP_sec >= 5 & doff_nossr(iD).x_LFP_sec <= 5.5;
            
            figure(1)
            clf
            subplot(2,1,1)
            plot(don_nossr(iD).x_LFP_sec(IX_ON), don_nossr(iD).LFP(IX_ON),'LineWidth',3)
            hold on
            plot(doff_nossr(iD).x_LFP_sec(IX_OFF), doff_nossr(iD).LFP(IX_OFF),'LineWidth',3)
            legend('WF On','WF Off'); legend boxoff
            title([ 'No SSR ' tit ])
            axis tight
            
            pubify_figure_axis
            %         ylabel('$\mu$V','Interpreter','Latex')
            ylabel('mV','Interpreter','Latex')
            xlabel('sec')
            
            
            IX_ON = don_ssr(iD).x_LFP_sec >= 5 & don_ssr(iD).x_LFP_sec <= 5.5;
            IX_ONlong = don_ssr(iD).x_LFP_sec >= 5 & don_ssr(iD).x_LFP_sec <= 10;
            
            IX_OFF = doff_ssr(iD).x_LFP_sec >= 5 & doff_ssr(iD).x_LFP_sec <= 5.5;
            
            subplot(2,1,2)
            plot(don_ssr(iD).x_LFP_sec(IX_ON), don_ssr(iD).LFP(IX_ON),'LineWidth',3)
            hold on
            plot(doff_ssr(iD).x_LFP_sec(IX_OFF), doff_ssr(iD).LFP(IX_OFF),'LineWidth',3)
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
            [OUT(iD).welchssr, OUT(iD).welchnossr, OUT(iD).welchpure] = DANA_plot_psd_ssr_on_off(don_ssr(iD).LFP(IX_ONlong), don_nossr(iD).LFP(IX_ONlong), META.fq_logspace, LFP_sFreq,Target_Fq,doff_nossr(iD).LFP(IX_ONlong));
            if saveit, saveas(gcf,fullfile(figdir,sprintf('welchssrnossr%s',tit)),savetype), end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if PLOT_WAVELET
                % plot spectral response
                IX_ON = don_nossr(iD).x_LFP_sec >= 5 & don_nossr(iD).x_LFP_sec <= 6;
                IX_OFF = doff_nossr(iD).x_LFP_sec >= 5 & doff_nossr(iD).x_LFP_sec <= 6;
                figure(2)
                clf
                subplot(2,1,1)
                DANA_plot_wavelet_and_scans( don_nossr(iD).LFP(IX_ON),don_nossr(iD).x_LFP_sec(IX_ON),...
                    LFP_sFreq,don_nossr(iD).EventRecsDown_sec,META.fq_logspace,wavelet_size_in_scan);
                title([tit ' No SSR, Scanning On'])
                subplot(2,1,2)
                DANA_plot_wavelet_and_scans( doff_nossr(iD).LFP(IX_OFF),doff_nossr(iD).x_LFP_sec(IX_OFF),...
                    LFP_sFreq,[],META.fq_logspace,wavelet_size_in_scan);
                title('No SSR, Scanning Off')
                if saveit, saveas(gcf,fullfile(figdir,sprintf('waveletNoSSR%s',tit)),savetype), end
            end
            % plot template match
            % but choose place so that it starts with offset.
            ix = find(don_nossr(iD).EventRecsDown_sec > 5,1,'first');
            st_sec = don_nossr(iD).EventRecsDown_sec(ix);
            IX_ON = don_nossr(iD).x_LFP_sec >= st_sec & don_nossr(iD).x_LFP_sec <= st_sec + 1;
            IX_OFF = doff_nossr(iD).x_LFP_sec >= st_sec & doff_nossr(iD).x_LFP_sec <= st_sec + 1;
            
            [LFP_target, INFO] = Artificial_LFP(LFP_sFreq, sum(IX_ON)/LFP_sFreq, Target_Fq, 0, 0 );
            fqs = [.4:.2:250];
            step_size_samples = LFP_sFreq/100;
            
%             [rmatch,rmatch_ideal, INFO]= SPEC_time_to_recover_signal(don_nossr(iD).LFP(IX_ON), LFP_sFreq, LFP_target,fqs,step_size_samples,true);
            
            [~,Pwv_on] = SPEC_waveletdecomp(META.fq_logspace,don_nossr(iD).LFP(IX_ON),LFP_sFreq, ...
                wavelet_size_in_scan);
            
            [~,Pwv_off] = SPEC_waveletdecomp(META.fq_logspace,doff_nossr(iD).LFP(IX_OFF),LFP_sFreq,...
                wavelet_size_in_scan);

            figure(10)
            clf
            [rmatch,x,mk,OUT(iD).wavelet_psd_nossr]=DANA_plot_template_match(Pwv_on,don_nossr(iD).x_LFP_sec(IX_ON),META.fq_logspace,LFP_sFreq, Pwv_off,don_nossr(iD).EventRecsDown_sec);
            [rmatch_pure,~,~,OUT(iD).wavelet_psd_pure]=DANA_plot_template_match(Pwv_off,doff_nossr(iD).x_LFP_sec(IX_OFF),META.fq_logspace,LFP_sFreq, Pwv_off,don_nossr(iD).EventRecsDown_sec);
            ix = find(rmatch_pure.^2 > .7,1,'first')
            ix/LFP_sFreq
            
            OUT(iD).wavelet_psd_xlabs = META.fq_logspace;
            title([tit ' Template match no SSR'])
            if saveit, saveas(gcf,fullfile(figdir,sprintf('template_match_nossr%s',tit)),savetype), end
            
            ix = find(rmatch_pure.^2 > .7,1,'first');
            fprintf('Detected %d Hz at %1.2f s \n', Target_Fq, ix/LFP_sFreq);
            
            OUT(iD).rmatch_no_ssr = rmatch;
            OUT(iD).rmatch_pure = rmatch_pure;
            OUT(iD).rmatch_x_sec_no_ssr= x;
            OUT(iD).rmatch_scan_times_s_no_ssr= mk;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if PLOT_WAVELET
                IX_ON = don_ssr(iD).x_LFP_sec >= 5 & don_ssr(iD).x_LFP_sec <= 6;
                IX_OFF = doff_ssr(iD).x_LFP_sec >= 5 & doff_ssr(iD).x_LFP_sec <= 6;
                
                figure(3)
                clf
                subplot(2,1,1)
                DANA_plot_wavelet_and_scans( don_ssr(iD).LFP(IX_ON),don_ssr(iD).x_LFP_sec(IX_ON),...
                    LFP_sFreq,don_ssr(iD).EventRecsDown_sec,META.fq_logspace,wavelet_size_in_scan);
                title([tit ' SSR, Scanning On'])
                subplot(2,1,2)
                DANA_plot_wavelet_and_scans( doff_ssr(iD).LFP(IX_OFF),doff_ssr(iD).x_LFP_sec(IX_OFF),...
                    LFP_sFreq,[],META.fq_logspace,wavelet_size_in_scan);
                title('SSR, Scanning Off')
                if saveit, saveas(gcf,fullfile(figdir,sprintf('waveletSSR%s',tit)),savetype), end
            end
            % plot template match
            IX_ON = don_ssr(iD).x_LFP_sec >= 5 & don_ssr(iD).x_LFP_sec <= 6;
            IX_OFF = doff_ssr(iD).x_LFP_sec >= 5 & doff_ssr(iD).x_LFP_sec <= 6;
            [~,Pwv_on] = SPEC_waveletdecomp(META.fq_logspace,don_ssr(iD).LFP(IX_ON),LFP_sFreq,wavelet_size_in_scan);
            [~,Pwv_off] = SPEC_waveletdecomp(META.fq_logspace,doff_ssr(iD).LFP(IX_OFF),LFP_sFreq,wavelet_size_in_scan);
            
            figure(11)
            clf
            [rmatch,x,mk,OUT(iD).wavelet_psd_ssr]=DANA_plot_template_match(Pwv_on,don_ssr(iD).x_LFP_sec(IX_ON),META.fq_logspace,LFP_sFreq, Pwv_off,don_ssr(iD).EventRecsDown_sec);
            title([tit ' Template match SSR'])
            if saveit, saveas(gcf,fullfile(figdir,sprintf('template_match_ssr_HF%s',tit)),savetype), end
            drawnow
            OUT(iD).rmatch_ssr = rmatch;
            OUT(iD).rmatch_x_sec_ssr= x;
            OUT(iD).rmatch_scan_times_s_ssr= mk;
            
            
            figure(12)
            FQ_IX = META.fq_logspace >= 600;
            clf
            [ OUT(iD).rmatch_ssr_high_fq ,OUT(iD).rmatch_ssr_high_fq_x_sec,OUT(iD).rmatch_ssr_high_scan_times_s,OUT(iD).wavelet_psd_ssr_HF]= ...
                DANA_plot_template_match(Pwv_on(FQ_IX,:),don_ssr(iD).x_LFP_sec(IX_ON),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_ssr(iD).EventRecsDown_sec);
            title([tit ' Template match SSR HF'])
            OUT(iD).rmatch_ssr_high_fq_HF_x = META.fq_logspace(FQ_IX);
            
            [OUT(iD).rmatch_nossr_HF ,OUT(iD).rmatch_nossr_HF_x_sec, OUT(iD).rmatch_nossr_HF_scan_times_s,OUT(iD).wavelet_psd_nossr_HF]= ...
                DANA_plot_template_match(Pwv_off(FQ_IX,:),don_nossr(iD).x_LFP_sec(IX_ON),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_ssr(iD).EventRecsDown_sec,false);
            
            [OUT(iD).rmatch_pure_HF ,OUT(iD).rmatch_pure_HF_x_sec, OUT(iD).rmatch_pure_HF_scan_times_s,OUT(iD).wavelet_psd_pure_HF]= ...
                DANA_plot_template_match(Pwv_off(FQ_IX,:),doff_nossr(iD).x_LFP_sec(IX_OFF),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_ssr(iD).EventRecsDown_sec,false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure(14)
            FQ_IX = META.fq_logspace < 250;
            clf
            [OUT(iD).rmatch_ssr_low_fq, OUT(iD).rmatch_ssr_low_fq_x_sec,mkLF,OUT(iD).wavelet_psd_ssr_LF]=DANA_plot_template_match(Pwv_on(FQ_IX,:),don_ssr(iD).x_LFP_sec(IX_ON),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_ssr(iD).EventRecsDown_sec);
            title([tit ' Template match SSR LF'])
            
            OUT(iD).rmatch_ssr_low_fq_HF_x = META.fq_logspace(FQ_IX);
            OUT(iD).rmatch_ssr_low_scan_times_s= mkLF;
            
            [OUT(iD).rmatch_nossr_LF ,OUT(iD).rmatch_nossr_LF_x_sec, OUT(iD).rmatch_nossr_LF_scan_times_s,OUT(iD).wavelet_psd_nossr_LF]= ...
                DANA_plot_template_match(Pwv_on(FQ_IX,:),don_nossr(iD).x_LFP_sec(IX_ON),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_ssr(iD).EventRecsDown_sec,false);
            
            [OUT(iD).rmatch_pure_LF ,OUT(iD).rmatch_pure_LF_x_sec, OUT(iD).rmatch_pure_LF_scan_times_s,OUT(iD).wavelet_psd_pure_LF]= ...
                DANA_plot_template_match(Pwv_off(FQ_IX,:),doff_nossr(iD).x_LFP_sec(IX_OFF),META.fq_logspace(FQ_IX),LFP_sFreq, Pwv_off(FQ_IX,:),don_nossr(iD).EventRecsDown_sec,false);
            
        OUT(iD).Target_Fq = Target_Fq;
        OUT(iD).LFP_sFreq = LFP_sFreq;
        
        fprintf('.')
    end
    % Save pre data.
    disp('Saving')
    [~,fn] = fileparts(fnames{iF});
    
    save([mname fn],'OUT','don_ssr','mname');
end
% load(mname)
%
%
%% Template matching SnR analysis.
savetype = 'eps';
FQ  = [don_ssr.Freq];
xIX = FQ >= 1;
LFP_sFreq = don_ssr(1).LFP_sFreq;
Rssr = [];
RssrHF = [];
Rnossr = [];
Rpure = [];
for ii = 1:length(OUT)
    Rssr(ii,:) =  OUT(ii).rmatch_ssr(1:9999)';
    Rnossr(ii,:) =  OUT(ii).rmatch_no_ssr(1:9999)';
    Rpure(ii,:) =  OUT(ii).rmatch_pure(1:9999)';
    RssrHF(ii,:) =  OUT(ii).rmatch_ssr_high_fq(1:9999)';
    RnossrHF(ii,:) =  OUT(ii).rmatch_nossr_HF(1:9999)';
    RpureHF(ii,:) =  OUT(ii).rmatch_pure_HF(1:9999)';
end
Rssr = Rssr(xIX,:);
Rnossr = Rnossr(xIX,:);
Rpure = Rpure(xIX,:);
RssrHF = RssrHF(xIX,:);
RnossrHF = RnossrHF(xIX,:);
RpureHF = RpureHF(xIX,:);
%%
newFQ = FQ(xIX);
ix = 1:LFP_sFreq/5;
flab = [];
for ii = 1:length(newFQ)
    flab{ii} = [num2str(newFQ(ii)) 'Hz'];
end


figure
subplot(1,3,1)
imagesc((1:length(ix))/LFP_sFreq ,[],Rnossr(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('No SSR')
set(gca,'XLim',[0 .18])


subplot(1,3,2)
imagesc((1:length(ix))/LFP_sFreq ,[],Rssr(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('SSR')
set(gca,'XLim',[0 .18])

subplot(1,3,3)
imagesc((1:length(ix))/LFP_sFreq ,[],Rpure(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')
title('Base: no Scans, No SSR ')
set(gca,'XLim',[0 .18])

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