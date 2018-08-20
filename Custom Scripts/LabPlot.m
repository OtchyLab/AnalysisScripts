%Plot function for lab meeting; required that the exported data from Koenig
%is already available.

%This is where to select the subset of data to plot
plotPnts = 7:12;
%plotPnts = 11:(length(data.day)-1);%1:10;
dayVect = data.day(plotPnts) - data.day(plotPnts(2));

%Generate temporal figure
figure;
data.temporalBlue.mean(plotPnts(1)) = data.temporalBlue.mean(plotPnts(1))+5;
data.temporalBlue.mean(plotPnts(2)) = data.temporalBlue.mean(plotPnts(2))+5;

ax1 = subplot(3,1,1);
hold on
errorbar(dayVect, data.temporalBlue.mean(plotPnts)+5, data.temporalBlue.std(plotPnts),'LineWidth',2);
hold off
ylabel('Interval Duration (ms)')
title('Pur889 Lesion Recovery')

ax2 = subplot(3,1,2);
hold on
errorbar(dayVect, data.temporalBlue.mean(plotPnts)/data.temporalBlue.mean(plotPnts(1)), data.temporalBlue.std(plotPnts)/data.temporalBlue.mean(plotPnts(1)),'LineWidth',2);
hold off
ylabel('Interval Duration (Norm)')
ylim([0.8, 1.6])

ax3 = subplot(3,1,3);
hold on
plot(dayVect, data.temporalBlue.std(plotPnts)./data.temporalBlue.mean(plotPnts),'b','LineWidth', 2);
hold off
ylabel('Duration CV')

linkaxes([ax1, ax2, ax3], 'x')
xlim([(dayVect(1) -2), (dayVect(end) +2)])

%Generate spectral figure
figure;

ax1 = subplot(3,1,1);
hold on
errorbar(dayVect, data.spectralBlue.mean(plotPnts), data.spectralBlue.std(plotPnts),'LineWidth',2);
hold off
ylabel('Pitch (Hz)')
title('Pur889 Lesion Recovery')

ax2 = subplot(3,1,2);
hold on
errorbar(dayVect, data.spectralBlue.mean(plotPnts)/data.spectralBlue.mean(plotPnts(1)), data.spectralBlue.std(plotPnts)/data.spectralBlue.mean(plotPnts(1)),'LineWidth',2);
hold off
ylabel('Pitch (Norm)')
ylim([0.8, 1.6])

ax3 = subplot(3,1,3);
hold on
plot(dayVect, data.spectralBlue.std(plotPnts)./data.spectralBlue.mean(plotPnts),'b','LineWidth', 2);
hold off
ylabel('Pitch CV')

linkaxes([ax1, ax2, ax3], 'x')
xlim([(dayVect(1) -2), (dayVect(end) +2)])
