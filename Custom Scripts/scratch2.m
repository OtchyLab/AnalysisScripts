%Plotting single shot TDT data
% pick = [2, 4, 5, 6, 7];
pick = 1:8;
d = washdata;

trials = 8;
chans = 6;

% plotColor = [0.5, 0.5, 0.5];
% plotColor = [1, 0, 0];
plotColor = [0, 1, 0];


 
%Single figure with all of the channels in subplots
figure(1001); %clf
figure(1002); %clf


for i = 2
   figure(1001); %subplot(3, 2, i)
   
   %find offset
   %delta = mean(mean(d.tdt.response(:,1:75,i)));
   delta = (mean(d.tdt.response(pick,1:75,i),2));
   
   %Plot the raw traces
   plot(d.tdt.times_aligned*1000, (d.tdt.response(pick,:,i)-delta)*10e3, 'Color', plotColor); hold on 
   
   %Format the figure
   axis tight
   xlim([-3, 4])
   ylim([-.000125, .000125]*10e3)
   set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', -2:1:2)
   ylabel('Evoked Response (mV)')
   xlabel('Time (ms)')
   
   figure(1002); %subplot(3, 2, i)
   
   %Plot the mean trace
   shadedErrorBar(d.tdt.times_aligned*1000, mean(d.tdt.response(pick,:,i)-delta)*10e3, std((d.tdt.response(pick,:,i)-delta)*10e3, 1, 1), 'g', 1); hold on 
   
   %Format the figure
   axis tight
   xlim([-3, 4])
   ylim([-.000125, .000125]*10e3)
   set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', -2:1:2)
   ylabel('Evoked Response (mV)')
   xlabel('Time (ms)')
end

figure(1001);
% title(gcf,'Raw response to 64uA stimulation')
set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 5.5, 4.5])

figure(1002); %clf
% title(gcf,'Mean +/- Std response to 64uA stimulation')
set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 5.5, 4.5])

