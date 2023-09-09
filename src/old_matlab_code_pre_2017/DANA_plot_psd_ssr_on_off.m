function [welchssr, welchnossr,welchpure]=DANA_plot_psd_ssr_on_off(LFPssr,LFPnossr,fqs,sFreq, target_fq, pureLFP)
if nargin < 6
    pureLFP = [];
end
logscale = true;
welchpure = [];
welch_win = 3*sFreq; % 3s window - get as much as you can. The bigger, the better the fq resolution. For this, it's OK
lw = 2;
[welchssr] = pwelch(LFPssr,welch_win,[],fqs,sFreq); % let welch decide.
[welchnossr] = pwelch(LFPnossr,welch_win,[],fqs,sFreq); % let welch decide.
ix = find(fqs >= target_fq,1,'first');
%%
cla
h = plot(fqs,welchssr,'b','LineWidth',lw);
hold on
plot(fqs,welchnossr,'k','LineWidth',lw);
if ~isempty(pureLFP)
    [welchpure] = pwelch(pureLFP,welch_win,[],fqs,sFreq); % let welch decide.
    plot(fqs,welchpure,'g','LineWidth',lw+2);
    plot(fqs,welchssr,'b','LineWidth',lw)
    plot(fqs,welchnossr,'k','LineWidth',lw);
end
legend('ssr','nossr','pure')
% set(gca,'XScale','log')
% uistack(h,'top')
axis tight
a = axis;
plot(fqs(ix),a(4)*1.1,'rv','MarkerSize',8,'MarkerFaceColor','r')
plot(fqs(ix),a(3),'r^','MarkerSize',8,'MarkerFaceColor','r')

if logscale
    set(gca,'YScale','log')
    ylabel('$$ log10 \mu V^2 $$','Interpreter','latex')
else
    ylabel('$$ \mu V^2 $$','Interpreter','latex')
end
if target_fq <  50
    xl = [0 80];
elseif  target_fq >=  50 && target_fq < 200
    xl = [20 300];
elseif  target_fq >=  200
    xl = [100 3000];
end
set(gca,'XLim',xl)
title([num2str(target_fq) ' Hz'])
box off
hold on
xlabel('Hz')
% pubify_figure_axis
% grid on

