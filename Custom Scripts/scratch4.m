%Plotting all trials of a current sweep on a single figure

% %Ask user which file to load; check for suffix
% [filename, pathname, filterindex] = uigetfile('C:\Users\Tim\Desktop\LLR20-2017-02-13\*.mat', 'Pick a set of file for sweep', 'MultiSelect', 'on');
% numFiles = numel(filename);
% 
% %Read in the data from each of the Gamry files
% % dataBlock = [];
% if numFiles > 1
%     %Sequentially read in the nultiple files
%     for i = 1:numFiles
%         fullName = [pathname, filesep, filename{i}];
%         load(fullName);
%         dataBlock(i) = data;
%     end
% else
%     %Read in a single file
%     fullName = [pathname, filesep, filename];
%     dataBlock = load(fullName);
% end

% Number of channels to process
chans = 6;
stimMax = 120;
stimMin = 0;

%Set figure
figure(1005); clf
set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 5.5, 4.5])
clr = jet(100);

for i = numFiles:-1:1
    
    for j = 1
        %Set location
        %subplot(chans, 1, j)
        caxis([stimMin stimMax])
        
        %find offset
        delta = (mean(dataBlock(i).tdt.response(:,1:75,j),2));
        
        %Find Stim intensity
        stim = dataBlock(i).stim.current_uA;
        a = interp1(linspace(stimMin,stimMax,100), clr(:,1), stim);
        b = interp1(linspace(stimMin,stimMax,100), clr(:,2), stim);
        c = interp1(linspace(stimMin,stimMax,100), clr(:,3), stim);
        col = [a, b, c];

        %Plot the raw traces
        plot(dataBlock(i).tdt.times_aligned*1000, (dataBlock(i).tdt.response(1:1:end,:,j)-delta(1:1:end))*10e3, 'Color', col); hold on
        
        %Format the figure
        axis tight
        xlim([-2, 4])
        ylim([-.00015, .0002].*10e3)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', -2:1:2)
        ylabel(['Evoked Response (mV)'])
        colormap(jet); colorbar('TickDirection', 'out', 'Ticks', 0:30:120)
    end
    
    
end

xlabel('Time (ms)')

