function TS_ephys_analysisSUM
% This script runs a per-bird anaysis on the TS MUA recordings from birds.
%The main input for this function is the output of the makeWaterfallHiRes.m script.

%Set the output file location

%Select the dataset to use
dataFolder = 'C:\Users\Tim\Desktop\Poster Datasets\';
dataName = 'LR81RY177 TS Recordings (Ch1).mat';
templateName = 'LR81RY177_170825_template345.mat';
% cond = 'hit';

%Load the data
load([dataFolder, templateName]);
load([dataFolder, dataName]);

%close all; %close the windows

h = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the windowed waterfall image

%Just in case the properly windowed matrix wasn't created
if ~isfield(data, 'neuropowAbs_winBreaks')

    winSize = 3; winStep = 1;
    neuropowAbs_winBreaks = [];
    t = [0, data.dayEnds];
    for i = 2:length(t)
        inxRend = (t(i-1)+1):t(i);
        [indx_mat] = windowTS(inxRend,winSize,winStep,'pad');
        winBreaks = [];
        for j = 1:size(indx_mat,1)
            trimIndx = indx_mat(j,~isnan(indx_mat(j, :)));
            winBreaks(j, :) = mean(data.neuropowAbs(trimIndx,:),1);
        end
        neuropowAbs_winBreaks = [neuropowAbs_winBreaks; winBreaks];
    end
    data.neuropowAbs_winBreaks = neuropowAbs_winBreaks;
    
    % Plot the Windowed Correlation
    [data.neuroCovWinBreaks, ~] = corrcoef(data.neuropowAbs_winBreaks');
end

%Generate Time vector from remaining files
timeVect = [];
for i = 1:length(data.files)
    timeVect(i) = getFileTime(data.files{i});
end
colors = [{'b'}, {'r'}, {'g'}, {'c'}, {'m'}, {'y'}];
h(1) = figure(1); clf;%template
imagesc(-template);
axis tight; axis xy; colormap(jet)
set(gca, 'Box', 'off',  'TickDir', 'out', 'LineWidth', 2, 'FontSize', 14)
set(gca, 'XTickLabels', [], 'YTick', [])

h(2) = figure(2); clf;
t = [0, data.dayEnds];
subplot(1,6,1:3); cla %HVC waterfall
imagesc(data.neuropowAbs_winBreaks);
% imagesc(data.neuropowAbs); 
hold on

%Add in the day markers
offset = -500;
color = [{'r'}, {'g'}, {'c'}, {'m'}, {'y'}];
for i = 1:length(data.dayEnds)
    plot([offset, offset], [t(i), t(i+1)], 'Color', colorTint(i, 0), 'LineWidth', 3, 'clip', 'off'); hold on
end

axis tight; axis ij; colormap(jet)
x = xlim; xlim([0,x(2)]);
xticks = 100:100:(x(2)*1000/44150);
xpos = round((xticks/1000)*44150);
ypos = 2000:2000:size(data.neuropowAbs_winBreaks,1);
set(gca, 'Box', 'off',  'TickDir', 'out', 'XTick', xpos, 'XTickLabels', xticks, 'YTick', ypos, 'LineWidth', 2, 'FontSize', 8)
xlabel('Time (ms)', 'FontSize', 8)
ylabel('Renditions', 'FontSize', 8)

subplot(1,6,4); cla; %running correlation
hold on
recovBlock = data.neuroCovWinBreaks(1:data.dayEnds(1), :);
recovBlock(recovBlock ==1) = NaN;
recovCorr = nanmean(recovBlock,1);
winSize = 51; winStep = 1;
recovCorrSnips = []; rendSnips = [];
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    recovCorrSnips{i-1} = nanmean(windowTS(recovCorr(rendSnips),winSize,winStep,'pad', 'boxcar')');
    plot(recovCorrSnips{i-1}, rendSnips, '-', 'LineWidth', 2, 'Color', colorTint(i-1, 0))
end
axis tight; axis ij;
xlim([0, 1]);
xlabel('Corr to Pre-Lesion', 'FontSize', 8)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0.25:0.25:1], 'YTick', [], 'LineWidth', 2, 'FontSize', 8)

subplot(1,6,5); cla; %mean HVC Activity
hold on
t = [0, data.dayEnds];
meanPow = mean(data.neuropowAbs_winBreaks, 2);
meanPowSnips = [];
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    meanPowSnips{i-1} = nanmean(windowTS(meanPow(rendSnips),winSize,winStep,'pad', 'boxcar')');
    plot(meanPowSnips{i-1}, rendSnips,  '-', 'LineWidth', 2, 'Color', colorTint(i-1, 0))
end
axis tight; axis ij;
%xlim([0, 0.1])
set(gca, 'Box', 'off',  'TickDir', 'out', 'YTick', [], 'LineWidth', 2, 'FontSize', 8)
xlabel('Mean HVC Activity (V^2)', 'FontSize', 8)

subplot(1,6,6); cla; %bout break length
hold on
minBreak = 100;
motifDiff = [diff(timeVect) * (24*60*60), 0]; %in seconds
shortIdx = motifDiff < minBreak; %less than 500ms break and its the same bout
longIdx = motifDiff > (10*60*60); %more than 10 hour break and it's nighttime
motifDiff(shortIdx | longIdx) = NaN;
durDiffSnips = [];
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    durDiffSnips{i-1} = log10(motifDiff(rendSnips));
%     plot(durDiffSnips{i-1}, rendSnips,  '.', 'LineWidth', 2, 'MarkerSize', 4, 'Color', colors{i})
    for j = 1:length(durDiffSnips{i-1})
        if ~isnan(durDiffSnips{i-1}(j))
            plot([log10(minBreak), durDiffSnips{i-1}(j)], [rendSnips(j), rendSnips(j)],  '-', 'LineWidth', 1,'Color', colorTint(i-1, 0))
        end
    end
end
axis tight; axis ij;
xlim([log10(minBreak), 4.5])
set(gca, 'Box', 'off',  'TickDir', 'out', 'XTick', 2:2:4, 'YTick', [], 'LineWidth', 2, 'FontSize', 8)
xlabel('Bout Break Duration (log(s))', 'FontSize', 8)

set(gcf, 'Units', 'Inches', 'Position', [0 0 5 6])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the mean power during renditions
h(3) = figure(3); clf
cla; hold on
t = [0, data.dayEnds];
meanPow = mean(data.neuropowAbs, 2);
meanPowSnips = [];
% meanPowSnips{1} = meanPow(1:data.dayEnds(1));
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    meanPowSnips{i-1} = meanPow(rendSnips);
    scatter(rendSnips, meanPowSnips{i-1}, '.')
end
xlim([1, length(meanPow)])
xpos = 2000:2000:length(meanPow);
set(gca, 'Box', 'off',  'TickDir', 'out', 'XTick', xpos, 'LineWidth', 2, 'FontSize', 14)
xlabel('Renditions', 'FontSize', 18)
ylabel('Mean HVC Activity (V^2)', 'FontSize', 18)

set(gcf, 'Units', 'Inches', 'Position', [0 0 5 4])

h(4) = figure(4); clf
cla; hold on
t = [0, data.dayEnds];
sampleNum = 50;

%Collect the mean Activity
h(4) = figure(4); clf
cla; hold on
t = [1, data.dayEnds];
sampleSize = 50;
activSamp = [];
xs = (1:numel(data.neuropowAbs(1,:)))*1000/44150;
for i = 1:(numel(t)-1)
    rendSnips = t(i):(t(i)+sampleSize);
    activSamp(i,:) = smooth(mean(data.neuropowAbs(rendSnips,:), 1), 441, 'moving');
    plot(xs, activSamp(i,:), 'Color', colorTint(i, 0),'LineWidth', 2); hold on
end

axis tight;
xlabel('Time (ms)'); ylabel('Mean Activity (V^2)')
set(gca, 'Box', 'off',  'TickDir', 'out', 'LineWidth', 2, 'FontSize', 14)
%set(gca, 'XTickLabels', [], 'YTick', [])




% 
% %Plot the check in points
% sampSize = 50;
% meanPowStatM = mean(meanPowSnips{1});
% meanPowStatS = std(meanPowSnips{1});
% for i = 2:length(meanPowSnips)
%     %From the beginning of the snip
%     meanPowStatM(end+1) = mean(meanPowSnips{i}(1:sampSize));
%     meanPowStatS(end+1) = std(meanPowSnips{i}(1:sampSize));
%     
%     %From the end of the snip
%     meanPowStatM(end+1) = mean(meanPowSnips{i}(end-sampSize:end));
%     meanPowStatS(end+1) = std(meanPowSnips{i}(end-sampSize:end));
% end
% % x.D = [1, 4, 6, 8];
% % x.N = [3, 5, 7, 9];
% % x.L = 2;
% % x.all = 1:9;
% % 
% % h(4) = figure(4); clf
% % cla; hold on
% % errorbar(x.D, meanPowStatM(x.D), meanPowStatS(x.D), 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
% % errorbar(x.N, meanPowStatM(x.N), meanPowStatS(x.N), 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% % errorbar(x.L, meanPowStatM(x.L), meanPowStatS(x.L), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
% % 
% % xlim([0.5, 9.5])
% % ylabel('Mean HVC Activity (V^2)', 'FontSize', 18)
% % set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
% % set(gca, 'XTick', x.all, 'XTickLabels', [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}])
% % 
% % set(gcf, 'Units', 'Inches', 'Position', [0 0 5 4])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Neural Dynamics recovery
% subSize = 25;
% WinCorrM = mean(recovCorrSnips{1});
% WinCorrS = std(recovCorrSnips{1});
% for i = 2:length(recovCorrSnips)
%     %Beginning of current snip
%     WinCorrM(end+1) = mean(recovCorrSnips{i}(1:subSize));
%     WinCorrS(end+1) = std(recovCorrSnips{i}(1:subSize));
%     
%     %End of current snip
%     WinCorrM(end+1) = mean(recovCorrSnips{i}(end-subSize:end));
%     WinCorrS(end+1) = std(recovCorrSnips{i}(end-subSize:end));
% end
% 
% % h(5) = figure(8); clf
% % cla;
% % hold on
% % x = [];
% % x.D = [1, 4, 6, 8];
% % x.N = [3, 5, 7, 9];
% % x.L = 2;
% % x.all = 1:9;
% % 
% % errorbar(x.D, WinCorrM(x.D), WinCorrS(x.D), 'sb', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
% % errorbar(x.N, WinCorrM(x.N), WinCorrS(x.N), 'sk', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% % errorbar(x.L, WinCorrM(x.L), WinCorrS(x.L), 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
% % xlim([0.5, 9.5])
% % ylabel('HVC Dynamics Recovery', 'FontSize', 14)
% % set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 12)
% % set(gca, 'XTick', x.all, 'XTickLabels', [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}], 'YTick', [0.5, 1])
% % 
% % set(gcf, 'Units', 'Inches');
% % set(gcf, 'Position', [0 0 5 4])
% 
% %New Neural Dynamics recovery
% %Assemble the matrix for correlations
% origWin = 25;
% winOffset = ceil(origWin/2);
% neuralMat = [];
% for i = 1:(length(data.dayEnds)-1)
%     beforeEnd = data.neuropowAbs_winBreaks(data.dayEnds(i)-winOffset, :);
%     afterEnd = data.neuropowAbs_winBreaks(data.dayEnds(i)+winOffset, :);
%     neuralMat = [neuralMat; beforeEnd; afterEnd];
% end
% beforeEnd = data.neuropowAbs_winBreaks(data.dayEnds(end)-winOffset, :);
% neuralMat = [neuralMat; beforeEnd];
% 
% %Take correlation of the pre-lesion with each subsequent rendition
% neuroCorr = [];
% for i = 2:size(neuralMat,1)
%     neuroCorr(i) = corr(neuralMat(1,:)', neuralMat(i,:)');
% end
% 
% %Caluclate the correlation between non-pverlapping pre-lesion windows
% prepre =  data.neuropowAbs_winBreaks(data.dayEnds(1)-2*origWin, :);
% neuroCorr(1) = corr(neuralMat(1,:)', prepre');
% 
% % h(6) = figure(9); clf
% % plot(x.all, neuroCorr, 'sk', 'MarkerSize', 8)
% % xlim([0.5, 9.5])
% % ylim([0, 1])
% % ylabel('HVC Dynamics Recovery 2', 'FontSize', 14)
% % set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 12)
% % set(gca, 'XTick', x.all, 'XTickLabels', [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}], 'YTick', [0.5, 1])
% % 
% % set(gcf, 'Units', 'Inches', 'Position', [0 0 5 4])
% 
% sp = regexp(templateName, '_', 'split');
% 
% %Save data to file
% outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
% outputFile = [char(sp{1}) '_longRec 05192015'];
%Save all data
% save([outputLocation, outputFile, '_allData.mat']);
% 
% %Save figures
% savefig(h, [outputLocation, outputFile, '_figures.fig']);
% 
% %%%%%%%%%%%%%%%%%%%%%
% % Save to output location
% %%%%%%%%%%%%%%%%%%%%%
% output = [];
% output.name = char(sp{1});
% output.type = cond;
% output.dataName = dataName;
% output.xs = x;
% output.labels =  [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}];
% 
% output.meanPowStatM = meanPowStatM;
% output.meanPowStatS = meanPowStatS;
% output.WinCorrM = WinCorrM;
% output.WinCorrS = WinCorrS;
% output.neuroCorr = neuroCorr;
% 
% outPerm =  ['Long Rec Summary 05192015.mat'];
% m = exist([outputLocation, outPerm]);
% if m ==2 
%     %File already exists
%     clear('longRecStats');
%     load([outputLocation, outPerm], 'longRecStats')
%     longRecStats(end+1) = output;
% else
%     %No file yet created
%     longRecStats = output;
% end
% 
% %Save the updated data to file
% save([outputLocation, outPerm], 'longRecStats')
% 
% close all
% 
% display('done')

function cT = colorTint(fileNum, rendFrac)
%Based on the file snumber and rendition fraction, return a tint/shade of a
%pre-specified color

%Define colors (hg2 standards)
colorBase = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];

%We don't want to get too light (i.e., white), so scale with max = 0.5
fraction = rendFrac * 0.5;

%Generate tint by scaling color vector
cT = (1-colorBase(fileNum, :)) * fraction + colorBase(fileNum, :);


