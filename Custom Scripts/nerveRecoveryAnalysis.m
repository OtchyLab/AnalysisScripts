function nerveRecoveryAnalysis
%Scratch paper for labmeeting analysis

mother = 'V:\SongbirdData\Analysis';
% folder = 'llb59rblk132';
% folder = 'lr11ry176';
% folder = 'lr14rblk17';
% folder = 'LLY12';
% folder = 'LLY13';
% folder = 'LLY14';
% folder = 'LLY16';
% folder = 'LLR04';
% folder = 'LLY72';
% folder = 'LW58';
% folder = 'LW60';
 folder = 'LY80';
sourceDir = [mother, filesep, folder];

datafiles = dir([sourceDir, filesep, '*.mat']);
files = getStructField(datafiles, 'name');
indArray = 1:numel(files);

%%%%%%%%%%%%%%%%%%
%%%Spectrograms
%%%%%%%%%%%%%%%%%%
%Generate the mean spectrograms for each of the StretchEm datasets
spectro = [];
boolI = contains(files, 'dataset.mat');
ind = indArray(boolI);
for i = 1:numel(ind)
    %Load the file
    fullPath = [sourceDir, filesep, char(files(ind(i)))];
    load(fullPath, 'audio');
    
    %Take the mean of the aligned spectrograms
    spectro(i).pic = squeeze(mean(audio.aligned_audioCube, 1));
    
    %Note the record date of the data
    sp = regexp(files(ind(i)), '_', 'split');
    spectro(i).date = sp{1,1}(1,2);
    
    %Delete the loaded variable from the workspace
    clear('audio')
end
birdname = sp{1,1}(1,1);

%Plot specrtrograms as a single figure
labels = {'Pre', 'Post', 'Recovery'};
h(1) = figure(100); clf
for i = 1:numel(spectro)
    %Subplot select)
    subplot(3,1,i)
    imagesc(-spectro(i).pic)
    axis xy; axis tight
    colormap(jet)
    
    %Format
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTickLabels', '', 'LineWidth', 1)
    ylabel(labels(i))
    if i == 1
        title(birdname)
    end
    if i ~= numel(spectro)
        set(gca, 'XTickLabels', '')       
    end
end
%labels
xlabel('Time (ms)')
set(gcf, 'Units', 'inches', 'Position', [1, 1, 4, 6])


% %%%%%%%%%%%%%%%%%%
% %%%Temporal Structure
% %%%%%%%%%%%%%%%%%%
% %Load the file
% boolI = contains(files, 'temporal');
% ind = indArray(boolI);
% fullPath = [sourceDir, filesep, char(files(ind))];
% load(fullPath, 'dataOut');
% 
%Day Index
days = [-3:-1, 1:11]; %Days for plotting
plotSymbols = {'x', 's', '+', 'd'};
% numTraces = numel(dataOut.dataType);
% 
% %Plot a timeseries of syl/motif durations
% normM = [];
% h(2) = figure(101); clf
% for i = 1:(numTraces-1)
%     m = size(dataOut.int_m, 1);
%     
%     %Raw
%     subplot(2,1,1)
%     errorbar(days(1:m), dataOut.int_m(:,i), dataOut.int_std(:,i), [':' plotSymbols{i}], 'LineWidth', 1)
%     hold on
%     
%     xlim([-4, 9])
%     ylabel('Duration (ms)')
%     
%     %Normalized to baseline
%     subplot(2,1,2)
%     normM(:,i) = dataOut.int_m(:,i)./dataOut.int_m(3,i);
%     normS(:,i) = dataOut.int_std(:,i)./dataOut.int_m(3,i);
%     errorbar(days(1:m), normM(:,i), normS(:,i), [':' plotSymbols{i}], 'LineWidth', 1)
%     hold on
%     
% end
% %Plot the mean normalized
% % errorbar(days(1:m), mean(normM, 2), std(normM, 1, 2), '-ok', 'LineWidth', 1)
% shadedErrorBar(days(1:m), mean(normM, 2), std(normM, 1, 2), 'k', 1)
% 
% %Format
% xlim([-4, 9])
% ylabel('Norm Duration'); xlabel('Time (days)')
% set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1)
% 
% subplot(2,1,1)
% title(birdname)
% set(gca, 'Box', 'off', 'TickDir', 'out','XTickLabels', '', 'LineWidth', 1)
% 
% %Size
% set(gcf, 'Units', 'inches', 'Position', [1, 1, 4, 6])
% 
% %Zoom in
% h(3) = figure(102); clf
% for i = 1:(numTraces-1)
%     %Normalized pre and post
%     errorbar([1,2], [1, normM(4,i)], [0, normS(4,i)], ['-' plotSymbols{i} 'k'], 'LineWidth', 1)
%     hold on
%     
% end
% bar([1,2], [1, mean(normM(4,2:(numTraces-1)))], 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'LineWidth', 1, 'BarWidth', 0.5) 
% %Format
% xlim([0.5, 2.5]); ylim([0.75, 1.25])
% title(birdname); ylabel('Norm Duration')
% set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2], 'XTickLabels', {'Pre', 'Post'}, 'LineWidth', 1)
% 
% %Size
% set(gcf, 'Units', 'inches', 'Position', [1, 1, 3, 2])

%Delete the loaded file
clear('dataOut')

%%%%%%%%%%%%%%%%%%
%%%Acoustic recovery
%%%%%%%%%%%%%%%%%%
%Load the file
boolI = contains(files, 'spectral');
ind = indArray(boolI);
fullPath = [sourceDir, filesep, char(files(ind))];
load(fullPath, 'dataOut');

%Plot a timeseries of syl/motif durations
numTraces = numel(dataOut.dataType);
%normM = [];
h(4) = figure(103); clf
for i = 1:(numTraces)
    m = size(dataOut.sylM, 1);
    
    errorbar(days(1:m), dataOut.sylM(:,i), dataOut.sylS(:,i), [':' plotSymbols{i}], 'LineWidth', 1)
    hold on
   
end
%Plot the mean normalized
%errorbar(days(1:m), mean(dataOut.sylM, 2), std(dataOut.sylM, 1, 2), '-ok', 'LineWidth', 1.5)
meanSylM = mean(dataOut.sylM, 2); stdSylM = std(dataOut.sylM, 1, 2);
shadedErrorBar(days(1:m), meanSylM, stdSylM, 'k', 1)

%Format
title(birdname)
xlabel('Time (days)')
ylabel('Similarity to Baseline')
xlim([-4, 9]); ylim([0.5, 1.2])

set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1)

%Size
set(gcf, 'Units', 'inches', 'Position', [1, 1, 4, 3])

%Delete the loaded file
clear('dataOut')

%%%%%%%%%%%%%%%%%%
%%%Saving functions
%%%%%%%%%%%%%%%%%%
%Save all
save([sourceDir, filesep, 'allProcessData.mat'])

%Save figures to file
savefig(h(1), [sourceDir, filesep, 'spectrograms.fig'])
% savefig(h(2), [sourceDir, filesep, 'temporal.fig'])
% savefig(h(3), [sourceDir, filesep, 'temporalZoom.fig'])
savefig(h(4), [sourceDir, filesep, 'spectral.fig'])

saveas(h(1), [sourceDir, filesep, 'spectrograms.png'])
% saveas(h(2), [sourceDir, filesep, 'temporal.png'])
% saveas(h(3), [sourceDir, filesep, 'temporalZoom.png'])
saveas(h(4), [sourceDir, filesep, 'spectral.png'])

%Save to a common file for joint printing
jointName = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\Functional Recovery Figure\nerveBirdsData.mat';

output = [];
output.birdname = birdname;
output.manipulation = 'intact';
output.days = days(1:m);
% output.normMeanTemp = mean(normM, 2); clear

% output.normStdTemp = std(normM, 1, 2);
output.normMeanTemp = []; 
output.normStdTemp = [];
output.MeanSpect = meanSylM;
output.StdSpect = stdSylM;

x = exist(jointName);
if x ==2 
    %File already exists
    clear('summaryData');
    load(jointName, 'summaryData')
    summaryData(end+1) = output;
else
    %No file yet created
    summaryData = output;
end

% Save the updated data to file
save(jointName, 'summaryData')

