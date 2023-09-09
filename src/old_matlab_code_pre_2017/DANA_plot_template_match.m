function [rmatch,x,mk,template] = DANA_plot_template_match(SPEC,x,y_fq,LFP_sFreq, SPEC_base,mk, PLOT_IT)
if nargin < 7
    PLOT_IT = true;
end
% GIX = round(LFP_sFreq/4):Cols(SPEC)-round(LFP_sFreq/4);
midix = round(Cols(SPEC_base)/2);
template = SPEC_base(:,midix);
rmatch = zeros(Cols(SPEC),1);
for ii = 1:Cols(SPEC)
    tmp = corrcoef(SPEC(:,ii),template(:));
    rmatch(ii) = tmp(2);
end
mk = mk(mk > x(1) & mk < x(end));

if PLOT_IT
    subplot(2,1,1)
    plot(x,rmatch,'LineWidth',4)
    axis tight
    hold on
    set(gca,'YLim',[-1 1])
    a = axis;
    plot(mk,repmat(a(4),1,length(mk)),'gv','MarkerSize',10,'MarkerFaceColor','g')
    hold on
    plot(mk,repmat(a(3),1,length(mk)),'g^','MarkerSize',10,'MarkerFaceColor','g')
    ylabel('r')
    title(['Match to baseline'])
    xlabel('Seconds')
    
    subplot(2,1,2)
    plot(y_fq,template,'LineWidth',4)
    axis tight
    title('No Scan: the template');
    xlabel('Hz')
    ylabel('Power mV')
    
    set(gca,'XScale','log')
    grid on
    
    subplot(2,1,1)
end