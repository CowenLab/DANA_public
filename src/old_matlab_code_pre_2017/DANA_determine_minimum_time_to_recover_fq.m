function [rmatch] = DANA_determine_minimum_time_to_recover_fq()
%%
mfile = 'DANA_determine_minimum_time_to_recover_fq';
Frequency_range_to_examine_LOW = [    1 0 40;
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
    200 40 300];

Frequency_range_to_examine_HI = [    500 80 700;
    1000 600 1800;
    2000 1500 2400;
    2857 2400 3600; ];

META.fq_logspace = unique((logspace(log10(1),log10(3000),1650)));
META.fq_logspace_HI = unique(META.fq_logspace(META.fq_logspace >= 600));
META.fq_logspace_LO = unique(META.fq_logspace(META.fq_logspace <= 250));

% fqs = [.2:.2:250];
% fqs = META.fq_logspace_LO ;
% LFP_sFreq = 600;
% method = 'wavelet_cohen';
method = 'wavelet';
% method = 'pwelch';
% method = 'pwelch'; % does not work - still don't know precisely why -
% there is an important mystery as to why it does not work.
type = 'hi';
% type = 'lo';
noise = 0;
rmatch = []; INFO = [];
if strcmpi(type,'hi');
    FR = Frequency_range_to_examine_HI(:,1);
    LFP_sFreq = 10000;
    fqs = META.fq_logspace_HI;
    step_size_samples = round(LFP_sFreq/100);
elseif strcmpi(type,'lo')
    FR = Frequency_range_to_examine_LOW(:,1);
    LFP_sFreq = 1000;
    fqs = [.5:.5:250];
    %     fqs = [1:1:250];
    step_size_samples = round(LFP_sFreq/100);
end
RM = []; IN = [];
for ii = 1:length(FR)
    
    tgt_Hz = FR(ii);
    
    [LFP_target] = Artificial_LFP(LFP_sFreq, 3, tgt_Hz, 0, 0 );
    [LFP] = Artificial_LFP(LFP_sFreq, 0.1915, tgt_Hz, 0, noise );
    %     figure
    %     plot(LFP_target,'b')
    %     hold on
    %     plot(LFP,'r')
    %     clf
    [~,RM(ii,:),IN{ii}]= SPEC_time_to_recover_signal(LFP, LFP_sFreq, LFP_target,fqs,step_size_samples,method,false);
    %     drawnow
    %     pause
end
%
OUT = [];
for ii = 1:Rows(RM)
    OUT.target_frequency(ii) = FR(ii);
    OUT.time_to_recover_ideal_s(ii) = IN{ii}.time_to_recover_ideal_s;
    OUT.time_to_recover_s(ii) = IN{ii}.time_to_recover_s;
    OUT.best_match_r(ii) = (RM(ii,end));
end
OUT.target_frequency =  OUT.target_frequency(:);
OUT.time_to_recover_ideal_s =  OUT.time_to_recover_ideal_s(:);
OUT.time_to_recover_s =  OUT.time_to_recover_s(:);
OUT.best_match_r =  OUT.best_match_r(:);
%
figure
subplot(3,1,1)
plot(OUT.target_frequency,OUT.time_to_recover_ideal_s*1000,'b')
hold on
plot(OUT.target_frequency,OUT.time_to_recover_ideal_s*1000,'^','MarkerSize',20,'MarkerFaceColor','k')
pubify_figure_axis
title(['time_to_recover  ' type])

subplot(3,1,2)
bar(1:length(OUT.time_to_recover_ideal_s),OUT.time_to_recover_ideal_s*1000)
hold on
set(gca,'XTickLabel',OUT.target_frequency)
xlabel('Hz');
ylabel('ms')
title('time_to_recover')
pubify_figure_axis

subplot(3,1,3)
bar(1:length(OUT.time_to_recover_ideal_s),OUT.best_match_r.^2)
hold on
set(gca,'XTickLabel',OUT.target_frequency)
xlabel('Hz');
ylabel('R^2')
title([ method ' EV of original signal'])
pubify_figure_axis

label_figure(mfile)

% save as an excel sheet.
T = struct2table(OUT);
writetable(T,[method '_' type '.csv'])
%%
saveas(gcf,[method '_' type])
save([method '_' type])
% Find earliest time for each frequency.



%%
% do it again for high frequencies.
%%%%%%%%%%%%%%%%%%%%%%%

LFP_sFreq = 10000;
step_size_samples = round(LFP_sFreq/100);
for ii = 1:Rows(Frequency_range_to_examine_HI)
    
    tgt_Hz = Frequency_range_to_examine_HI(ii,1);
    
    [LFP_target] = Artificial_LFP(LFP_sFreq, 1, tgt_Hz, 0, 0 );
    [LFP] = Artificial_LFP(LFP_sFreq, 1, tgt_Hz, 0, noise );
    
    figure
    plot(LFP_target,'b')
    hold on
    plot(LFP,'r')
    figure
    [~,RMhf(ii,:),INhf{ii}]= SPEC_time_to_recover_signal(LFP, LFP_sFreq, LFP_target,fqs,step_size_samples,method,true);
    drawnow
    
end


for ii = 1:length(IN)
    tgt(ii) = Frequency_range_to_examine_HI(ii,1);
    ttr(ii) = IN{ii}.time_to_recover_ideal_s;
    ttr2(ii) = IN{ii}.time_to_recover_s;
end
figure
plot(tgt,ttr*1000,'b')
hold on
plot(tgt,ttr*1000,'^','MarkerSize',20,'MarkerFaceColor','k')
hold on 

title([method ' Minmal Time to Recover at .7 R^2'])
ylabel('ms')
xlabel('Hz')
pubify_figure_axis
saveas(gcf,[method 'HI'])
save([method 'HI'])


