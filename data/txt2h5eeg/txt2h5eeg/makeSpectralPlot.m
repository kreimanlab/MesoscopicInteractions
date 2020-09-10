function [ ] = makeSpectralPlot( t, V, V_before, e, SAMPLE_RATE, addText, directory)

%makeSpectralPlot plots time and power-frequency series

% Constants
%SAMPLE_RATE = 500; %512; % 1/sec
SEGMENT_LENGTH = 10*SAMPLE_RATE; % number of samples, window length in pwelch()
OVERLAP = 0.1*SEGMENT_LENGTH; % number of samples for overlap
NFFT = 1*SAMPLE_RATE; % number of DFT calculations

% Plot
h1 = figure;
%set(gcf,'Position',[0 300 600 750]);


% first subplot: time series of after filtering
subplot(3,1,2)
plot(t,V,'black')
title(strcat('Electrode',num2str(e),' time series after preprocessing'))
xlabel('Time (ms)')
ylabel('Voltage (uV)')
%axis([min(t) max(t) min(V) max(V)])
axis tight

% plot psd of after filtering
subplot(3,1,3)
[pxx,f] = pwelch(V(~isnan(V)),SEGMENT_LENGTH,OVERLAP,NFFT,SAMPLE_RATE);
pxxl1 = 10*log10(pxx);
plot(f,pxxl1,'black')
%semilogx(f,pxxl1,'black')
title(strcat('Electrode',num2str(e),' PSD after preprocessing'))
xlabel('Frequency (Hz)')
ylabel('Power/frequency (dB/Hz)')
%axis([min(f) max(f) min(pxxl1) max(pxxl1)])
%axis([1 max(f) min(pxxl1) max(pxxl1)])
axis tight
hold on
pmarker = [min(pxxl1) max(pxxl1)];
mfRange = 60:60:(SAMPLE_RATE/2);
mfLegend = cell(1,length(mfRange)+1);
mfLegend{1} = 'psd';
mfi = 2;
for mf = mfRange
    plot([mf mf],pmarker,'red--')
    mfLegend{mfi} = sprintf('%iHz',mf);
    mfi = mfi + 1;
end
legend(mfLegend,'Location','northoutside','Orientation','Horizontal')
%plot([60 60],pmarker,'red--')
%plot([120 120],pmarker,'magenta--')
%plot([180 180],pmarker,'blue--')
%plot([240 240],pmarker,'green--')
%plot([212 212],pmarker,'yellow--')
%plot([152 152],pmarker,'cyan--')
%plot([1/1.06 1/1.06],pmarker,'black--') %why was this here again? Heart Rate?
%legend('psd','60Hz','120Hz','180Hz','240Hz','300Hz','360Hz',...
%    'Location','northoutside','Orientation','Horizontal')

% plot psd before any processing
subplot(3,1,1)
[pxx,f] = pwelch(V_before(~isnan(V_before)),SEGMENT_LENGTH,OVERLAP,NFFT,SAMPLE_RATE);
pxxl1 = 10*log10(pxx);
plot(f,pxxl1,'black')
%semilogx(f,pxxl1,'black')
title(strcat('Electrode',num2str(e),' PSD before preprocessing'))
xlabel('Frequency (Hz)')
ylabel('Power/frequency (dB/Hz)')
%axis([min(f) max(f) min(pxxl1) max(pxxl1)])
%axis([1 max(f) min(pxxl1) max(pxxl1)])
axis tight
hold on
pmarker = [min(pxxl1) max(pxxl1)];
mfRange = 60:60:(SAMPLE_RATE/2);
mfLegend = cell(1,length(mfRange)+1);
mfLegend{1} = 'psd';
mfi = 2;
for mf = mfRange
    plot([mf mf],pmarker,'red--')
    mfLegend{mfi} = sprintf('%iHz',mf);
    mfi = mfi + 1;
end
legend(mfLegend,'Location','northoutside','Orientation','Horizontal')


%plot([60 60],pmarker,'red--')
%plot([120 120],pmarker,'magenta--')
%plot([180 180],pmarker,'blue--')
%plot([240 240],pmarker,'green--')
%plot([212 212],pmarker,'yellow--')
%plot([152 152],pmarker,'cyan--')
%plot([1/1.06 1/1.06],pmarker,'black--') %why was this here again? Heart Rate?
%legend('psd','60Hz','120Hz','180Hz','240Hz','300Hz','360Hz',...
%    'Location','northoutside','Orientation','Horizontal')

% Save the figure only if addText is not null
if ~strcmp(addText,'null')
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 15];
    print(h1,strcat([directory '/electrode'],num2str(e),'_psd',addText),'-dpng','-r0')
end

end

