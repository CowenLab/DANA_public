function DANA_analyze_brain_gel_ANA(fname)
% load('DANA_analyze_brain_gel.mat')
% mname = 'DANA_analyze_brain_gel_ANA';
load(fname)
set(0, 'defaultAxesFontName', 'Arial');

mname = 'DANA_analyze_brain_gel_ANA';
%%
FQ  = [OUT.Target_Fq];
xIX = FQ >= 1;

%% Template matching SnR analysis.
R = [];
for ii = 1:length(OUT)
      R(ii,:) =  OUT(ii).rmatch_wv_mn;
      %      R(ii,:) =  OUT(ii).rmatch_mn;
end
R = R(xIX,:);
newFQ = FQ(xIX);
OUT(1).LFP_sFreq = 10000;
ix = 1:OUT(1).LFP_sFreq/5;

figure
subplot(1,2,1)
imagesc((1:length(ix))/OUT(1).LFP_sFreq ,[],R(:,ix).^2)
set(gca,'YTick',1:length(newFQ))
set(gca,'YTickLabel',newFQ)
h = colorbar();
h.Label.String = 'r^2';
colormap(jet);
xlabel('Seconds')
ylabel('Hz')

subplot(1,2,2)
set(0,'DefaultAxesColorOrder',magma(length(newFQ)))

flab = [];
for ii = 1:length(newFQ)
    flab{ii} = [num2str(newFQ(ii)) 'Hz'];
end
plot((1:length(ix))/OUT(1).LFP_sFreq,R(:,ix)'.^2,'Linewidth',2)
label_figure(mname);
label_figure(fname,'bottom left');
legend(flab)
legend boxoff
xlabel('Seconds')
ylabel('r^2')

saveas(gcf,'FtemplateSnR.eps')
%%
SnR_no_FSCV = [OUT.SnR_before_welch];
SnR_no_FSCV = SnR_no_FSCV(xIX);
SnR = [OUT.Mean_Welch_SnR_after];
SnR = SnR(xIX);
newFQ = FQ(xIX);

good_ix = find(SnR >=1.5);
y = log(SnR); lab = 'log(SnR)';

%%% PLOT SUMMARY OF SnR %%%
clr = lines(2);
figure
clf
bar(1:length(y),y,'k')
axis tight
hold on
hold on
tmp = SnR_no_FSCV;
tmp(end+1) = tmp(end);
stairs((1:length(tmp))-.5,log(tmp),'Color',clr(2,:),'LineWidth',3)
% plot(good_ix, zeros(size(good_ix)),'g^')
set(gca,'XTick',1:length(SnR))
set(gca,'XTickLabel',newFQ)
set(gca,'XTickLabelRotation',45)
% plot_horiz_line_at_zero();
plot_horiz_line_at_zero(log(1.5),[],[.5 .5 .5]);
hold on
% plot_horiz_line_at_zero(log(1.5));
legend('Snr','Snr no FSCV')
text(length(SnR)+.4,log(1.5),'<- 1.5')
pubify_figure_axis
ylabel(lab)
xlabel('Frequency (Hz)')
grid on
label_figure(mname);
saveas(gcf,'F1.eps')
%% % PLOT Waterfall of power spectra %%%

fqs_to_use = [1, 5 , 10, 50, 100, 1000];
FQ = [OUT.Target_Fq];
[f,fqix] = intersect(FQ,fqs_to_use);
nc = ceil(length(fqs_to_use)/2);
h = [];
figure
clf
for ii = 1:length(fqs_to_use)
    fq = fqs_to_use(ii);
    fq_ix = fqix(ii);
    
    h(ii) = subplot(2,nc,ii);
    plot( OUT(fq_ix).welch_fq, (OUT(fq_ix).Mean_Full_Welch_after)','k');
    axis tight
    if ii == 1
        %        ylabel('10log10 power')
        ylabel('power (uV^2)')
    end
    title([num2str(OUT(fq_ix).Target_Fq) ' Hz'])
    ix = find(OUT(fq_ix).welch_fq >= OUT(fq_ix).Target_Fq,1,'first');
    hold on
    plot(OUT(fq_ix).welch_fq(ix),(OUT(fq_ix).Mean_Full_Welch_after(ix)*1.1),'rv','MarkerSize',4)
    %    pubify_figure_axis
    box off
    hold on
    xlabel('Hz')
end
equalize_y_axes(h)
label_figure(mname);
saveas(gcf,'./Figures/F2.eps')

%% Time to reach a SNR threshold
SNR = [OUT.SNR_by_time]';
SNR = SNR(xIX,:)
x_ms = OUT(2).SNR_by_time_x_sec*1000;

figure
clf
imagesc(x_ms,1:Rows(SNR), log(SNR))
hold on
firsts = [];
for ii = 1:Rows(SNR)
    ix = find(SNR(ii,:) > 2);
    firsts(ii,1) = ii;
    if ~isempty(ix)
        plot(x_ms(ix),ii,'w+')
        firsts(ii,2) = x_ms(ix(1));
    else
        firsts(ii,2) = nan;
    end
end
firsts = firsts(~isnan(firsts(:,2)),:);
plot(firsts(:,2), firsts(:,1),'w','LineWidth',4)

colormap(jet)

xlabel('ms');ylabel('Frequency (Hz)')
set(gca,'YTick',1:Rows(SNR))
set(gca,'YTickLabel',newFQ)
colorbar_label('Log SnR')
pubify_figure_axis
label_figure(mname);
saveas(gcf,'F3.eps')

%% Same data but simpler...
figure
clf
plot(firsts(:,2), firsts(:,1),'k.-','LineWidth',4,'MarkerSize',44)
axis tight
a = axis;
a(3) = 1;
a(2) = 700;
a(1) = 0;
axis(a)
set(gca,'YTick',1:Rows(SNR))
set(gca,'YTickLabel',newFQ)
axis ij
xlabel('ms');ylabel('Frequency (Hz)')
pubify_figure_axis
label_figure([mname ' ' F_time_to_id_frequency]);
grid on

saveas(gcf,'F_time_to_id_frequency.eps')

%% Attenuation as a function of frequency. 
ybefore = [OUT.Peak_mV_from_hilbert_before];
yafter  = [OUT.Peak_mV_from_hilbert_after];
ybefore = ybefore(xIX);
yafter = yafter(xIX);

figure
clf
h = [];
h(1) = subplot(1,4,1:2);
bar(1:length(ybefore),ybefore,'k')
axis tight
hold on
% plot(good_ix, zeros(size(good_ix)),'g^')
set(gca,'XTick',1:length(ybefore))
set(gca,'XTickLabel',newFQ)
set(gca,'XTickLabelRotation',45)

% rotateticklabel(gca,45)

ylabel('Peak Voltage (uV)')
xlabel('Frequency (Hz)')
title('Before FSCV')
% pubify_figure_axis
box on
% plot_horiz_line_at_zero();
plot_horiz_line_at_zero(log(1.5),[],[.5 .5 .5]);
hold on
grid on

h(2) = subplot(1,4,3:4);
bar(1:length(yafter),yafter,'k')
axis tight
hold on
% plot(good_ix, zeros(size(good_ix)),'g^')
set(gca,'XTick',1:length(yafter))
set(gca,'XTickLabel',newFQ)
set(gca,'XTickLabelRotation',45)

ylabel('Peak Voltage (uV)')
xlabel('Frequency (Hz)')
% plot_horiz_line_at_zero();
plot_horiz_line_at_zero(log(1.5),[],[.5 .5 .5]);
% pubify_figure_axis
box on
hold on
grid on
title('after FSCV')
label_figure(mname);
equalize_y_axes(h);
saveas(gcf,'F5.eps')

%% As above, but compute the post-ssr relative to the pre-SSR levels so that
% pre-ssr attentuation is accounted for.
figure(6)
clf
bar(1:length(ybefore),yafter./ybefore,'k')
axis tight
hold on
% plot(good_ix, zeros(size(good_ix)),'g^')
set(gca,'XTick',1:length(ybefore))
set(gca,'XTickLabel',newFQ)
set(gca,'XTickLabelRotation',45)

xlabel('Frequency (Hz)')
ylabel('FrequencySSR/FrequencyNoSSR')
title('Attentuation relative to pre-SSR')
% plot_horiz_line_at_zero();
plot_horiz_line_at_zero(1);
pubify_figure_axis
hold on
grid on
label_figure(mname);
saveas(gcf,'F6.eps')
%% Create a nice function of this.
figure(7)
clf
rng = newFQ(1):2:newFQ(end);
newy = interp1(newFQ,yafter./ybefore,rng,'spline');
plot(rng,newy,'k','LineWidth',5)
axis tight
hold on
% plot(good_ix, zeros(size(good_ix)),'g^')
xlabel('Frequency (Hz)')
% set(gca,'XTickLabelRotation',45)

ylabel('FrequencySSR/FrequencyNoSSR')
title('Attentuation relative to pre-SSR')
% plot_horiz_line_at_zero();
plot_horiz_line_at_zero(1);
hold on
grid on
pubify_figure_axis
label_figure(mname);
saveas(gcf,'F7.eps')
