% function DANA_time_to_recover_across_all_electrodes()
%% DANA - time to recover - averaged across electrodes...
% aggregate relevant data..
d = dir('DANA_analyze_brain_gel_post_4_14_16_v2aamp-D-*');
HF_ONLY = false; fq_ix = 1:14;
% HF_ONLY = true; fq_ix = 15:17;
mfile = 'DANA_time_to_recover_across_all_electrodes';
all_all_first = [];
for iF = 1:length(d)
    load (d(iF).name);
    LFP_sFreq = don_ssr(1).LFP_sFreq;
    FQ  = [don_ssr.Freq];

    fIX = FQ >= 1; 
    
    newFQ = FQ(fIX);
    ix = 1:LFP_sFreq/5;
    flab = [];
    for ii = 1:length(newFQ)
        flab{ii} = [num2str(newFQ(ii)) 'Hz'];
    end
    
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
    Rssr = Rssr(fIX,:);
    Rnossr = Rnossr(fIX,:);
    Rpure = Rpure(fIX,:);
    RssrHF = RssrHF(fIX,:);
    RnossrHF = RnossrHF(fIX,:);
    RpureHF = RpureHF(fIX,:);
    %
    %
    % Compute the onset times...
    %
    %
    R = [];
    % fq_ix = 1:17;
    flabs = flab(fq_ix);
    R{1} = Rssr(fq_ix,:).^2;
    R{2} = Rnossr(fq_ix,:).^2;
    R{3} = Rpure(fq_ix,:).^2;
    %   x_s = OUT(1).rmatch_x_sec_ssr(1:end-1);
    %  x_s = OUT(1).rmatch_ssr_high_fq_x_sec(1:end-1);
    if HF_ONLY
        % Special case: Just look at High-frequencies and time to recover.
        R = [];
        R{1} = RssrHF(fq_ix,:).^2;
        R{2} = RnossrHF(fq_ix,:).^2;
        R{3} = RpureHF(fq_ix,:).^2;
    end
    x_s = OUT(1).rmatch_ssr_high_fq_x_sec(1:end-1);
    %%
    x_s = x_s - x_s(1);
    th = 0.7;
    clrs = hsv(length(R));
    axis_meth = 2;
    all_first = [];
    
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
        
        all_first(:,iR) = firsts(:,2);
        
    end
    all_all_first(:,:,iF) = all_first;
end
%
all_all_first = all_all_first(:,:,[1:4 7 8]);

mn_all_firsts = nanmean(all_all_first,3);
se_all_firsts = Sem(all_all_first,3);

if HF_ONLY
    csvwrite('time_to_recover_HF.csv',mn_all_firsts)
    csvwrite('time_to_recover_HFse.csv',se_all_firsts)
    csvwrite('time_to_recover_HF_raw.csv',all_first)

else
    csvwrite('time_to_recover_LF.csv',mn_all_firsts)
    csvwrite('time_to_recover_LFse.csv',se_all_firsts)
    csvwrite('time_to_recover_LF_raw.csv',all_first)
end


%% Now we need the time to recover plot.

figure(10101)
clf
% axis tight
for iR = 1:Cols(mn_all_firsts)
    plot(1:Rows(mn_all_firsts), mn_all_firsts(:,iR),'.-','LineWidth',4,'MarkerSize',25,'Color',clrs(iR,:))
    hold on
    errorb(1:Rows(mn_all_firsts),  mn_all_firsts(:,iR), se_all_firsts(:,iR))
end
hold on
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
%
label_figure(mfile)

figure
imagesc(squeeze(all_all_first(:,1,:)))

all_all_firsts = all_all_first(:,1,[1:4 7 8]);
