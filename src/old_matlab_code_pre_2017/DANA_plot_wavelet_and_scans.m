function Pwv_on2=DANA_plot_wavelet_and_scans(L,x_sec,sFreq,marker_times_s,wv_fqs,wavelet_size_in_scan)
% ASSUMES L is in mV
[~,Pwv_on2] = SPEC_waveletdecomp(wv_fqs,L,sFreq,wavelet_size_in_scan);
x_sec = x_sec - x_sec(1);

imagesc(x_sec,[],Pwv_on2)
set(gca,'YTick',1:6:Cols(Pwv_on2))
set(gca,'YTickLabel',wv_fqs(1:6:end))
colormap (jet)
axis xy
colorbar
% colorbar_label('dB')
ylabel('Hz')
hold on
if ~isempty(marker_times_s)
    marker_times_s = marker_times_s - x_sec(1);
    a = axis;
    plot(marker_times_s,repmat(a(4),1,length(marker_times_s)),'gv','MarkerSize',10,'MarkerFaceColor','g')
    hold on
    plot(marker_times_s,repmat(a(3),1,length(marker_times_s)),'g^','MarkerSize',10,'MarkerFaceColor','g')
%     plot_markers_simple(marker_times_s)
end
yyaxis right
plot(x_sec,L,'y')
set(gca,'YLim',[-.8 5.8])
xlabel('Seconds')
ylabel('mV')
set(gca,'FontSize',8)
