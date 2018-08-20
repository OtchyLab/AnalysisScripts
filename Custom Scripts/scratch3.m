%Plotting sepoarately all trials of a single sTDT recording

d = spondata;

trials = 5;
chans = 6;

% plotColor = [0.5, 0.5, 0.5];
% plotColor = [1, 0, 0];
% plotColor = [1, 1, 1];

%Single figure with all of the channels in subplots
figure(1003); clf
set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 4.5, 9.5])
% figure(1002); %clf


for i = 1:trials
   figure(1003); subplot(trials, 1, i)
   
   for j = 1:chans
        %find offset
        delta = mean(mean(d.tdt.response(:,103:201,j)));
   
        %Plot the raw traces
        plot(d.tdt.times_aligned*1000, d.tdt.response(i,:,j)-delta+((j-1)*1e-4)); hold on 
   end
   
   %Format the figure
   if i ~= trials
       set(gca, 'XTickLabels', '');
   end
   axis tight
   ylim([-1.5e-4, .0007])
   set(gca, 'Box', 'off', 'TickDir', 'out')
   ylabel(['Trial ' num2str(i), ' (V)'])
   
end
xlabel(['Time (ms)'])

subplot(trials, 1, 1); title('Spontaneous TS Activity')





% figure(1001);
% % title(gcf,'Raw response to 64uA stimulation')
% 
% 
% figure(1002); %clf
% % title(gcf,'Mean +/- Std response to 64uA stimulation')
% set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 6, 9])

