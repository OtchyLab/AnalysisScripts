function TS_ephys_analysis
% The script creates a waterfall plot of MUA recordings that have been aligned with StretchEM.
% It solely looks at the neurtal recording data... use LongSong to look at the song data for this experiment.
%
%Original draft: 09/11/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Source Datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceDir = 'C:\Users\Tim\Desktop\Poster Datasets\';
annotFilelist = {'LR81RY177_170825_dataset345R_ch';...
                        'LR81RY177_170830_dataset345R_ch';...
                        'LR81RY177_170904_dataset345R_ch';...
                        'LR81RY177_170907_dataset345R_ch'};
chan = '1';

%Filename based processing
sp = regexp(annotFilelist{1}, '_', 'Split');
birdname = sp{1};
figTitle = [birdname ' TS Recordings (Ch' chan ')'];


% Predefine variables
files = []; neuropow = []; neuropowAbs = []; dayEnds = []; dates = [];
                    
% Cycle through the files, loading into memory what is needed and deleting the rest
numAnnots = numel(annotFilelist);
trimLength = 10; %in ms
trim = floor((trimLength/1000)*44150);
for i = 1:numAnnots
    % Clean load from drive
    filenames = []; data = []; neuro = [];
    fname = [annotFilelist{i}, chan, '.mat'];
    load([sourceDir, fname], 'filenames', 'data', 'neuro');
    
    %Filter out the files flaged during StretchEm inspection
    flags = ~data.flag;
    filenames = filenames(flags);
%     neuro.aligned_abs = neuro.aligned_abs(flags, :);
    neuro.aligned_abs = neuro.aligned(flags, :);
    
    %Filter out noisy recordings
    meanrawPow = mean(neuro.aligned_abs, 2);

    %shortwave
    if numel(meanrawPow) > 100
        mean200 = nanmean(windowTS(meanrawPow, 201, 1, 'pad','boxcar')');
        std200 = nanstd(windowTS(meanrawPow, 201, 1, 'pad','boxcar')');
        threshH = mean200 + 0.5*std200;
        threshL = mean200 - 0.5*std200;
    else
        mean200 = nanmean(windowTS(meanrawPow, 51, 1, 'pad','boxcar')');
        std200 = nanstd(windowTS(meanrawPow, 51, 1, 'pad','boxcar')');
        threshH = mean200 + 0.5*std200;
        threshL = mean200 - 0.5*std200;
    end
    
    cleanIdx = (meanrawPow < threshH') & (meanrawPow > threshL');
    
    %Add relevant data to the running stacks
    files = [files, filenames(cleanIdx)];
    neuropowAbs = [neuropowAbs; neuro.aligned_abs(cleanIdx, trim:(end-trim))];
    
    %Get date string
    sp = regexp(annotFilelist{i}, '_', 'split');
    dates = [dates; sp(2)];
end

%Longwave
meantrimPow = mean(neuropowAbs,2);
mean1000 = nanmean(windowTS(meantrimPow, 1001,1,'pad','boxcar')');
std1000 = nanstd(windowTS(meantrimPow, 1001,1,'pad','boxcar')');
thresh2 = mean1000 + 2*std1000;
cleanIdxLong = meantrimPow<thresh2';

%apply
files = files(cleanIdxLong);
neuropowAbs = neuropowAbs(cleanIdxLong, :);

%Generate Time vector from remaining files
timeVect = [];
for i = 1:length(files)
    timeVect(i) = getFileTime(files{i});
end

%Find indices of the day/nights and lesion time.
ind = 1:length(files);
for i = 1:length(dates)
    annotNum = datenum([dates{i}, '235959'], 'yymmddHHMMSS');
    tmp =ind(timeVect<annotNum);
    dayEnds(i) = tmp(end);
end

% Add lesion marker
% tmp = ind(timeVect<critTime);
% lesionNum = tmp(end);
% dayEnds = sort([dayEnds, lesionNum]);

%Calculate the mean power
meanPow = mean(neuropowAbs, 2);

% Generate covariance matrices using windowing
% winSize = 25; winStep = 1;
winSize = 11; winStep = 1;
neuropowAbs_winBreaks = [];
t = unique(([0, dayEnds]));
for i = 2:length(t)
    inxRend = (t(i-1)+1):t(i);
    [indx_mat] = windowTS(inxRend,winSize,winStep,'pad');
    winBreaks = [];
    for j = 1:size(indx_mat,1)
        trimIndx = indx_mat(j,~isnan(indx_mat(j, :)));
        winBreaks(j, :) = mean(neuropowAbs(trimIndx,:),1);
    end
    neuropowAbs_winBreaks = [neuropowAbs_winBreaks; winBreaks];
end

[neuroCovWinBreaks, ~] = corrcoef(neuropowAbs_winBreaks');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = [];

%Plot the smoothed waterfall image
h(1) = figure(1); clf;
subplot(1,5,1:3); cla;
imagesc(neuropowAbs_winBreaks); % waterfall
hold on

%Format the subplot
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 12)
xlabel('Time (ms)', 'FontSize', 12)
ylabel('Rendition', 'FontSize', 12)
title(figTitle, 'FontSize', 12);

subplot(1,5,4); cla;
% Day markers
x = xlim;
t = [0, dayEnds];
color = {'b', 'r', 'g', 'k', 'c', 'm', 'y', 'b', 'r', 'g', 'k'};
for i = 1:length(dayEnds)
    patch([0, 0, 1, 1], [t(i), t(i+1), t(i+1), t(i)], color{i}); hold on
end
axis tight; axis ij
set(gca, 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [0 0 5 11]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the mean power during renditions
h(2) = figure(2);
cla;
hold on
t = [0, dayEnds];
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    scatter(rendSnips, meanPow(rendSnips), '.')
end

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Mean HVC Activity (V^2)', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the windowed correlation
h(3) = figure(3);
cla;

%Show data on axes
imagesc(neuroCovWinBreaks, [0, 1]);
hold on

% Day markers
x = xlim;
y = ylim;
t = [0, dayEnds];
for i = 2:length(t)
    plot([x(1) - .03*x(2), x(1) - .03*x(2)], [(t(i-1)+8), t(i)-8],'LineWidth', 4, 'clip', 'off')
end

%Format the image
box off
set(gca, 'XTick', [], 'YTick', [])
xlabel('Increasing Time \Rightarrow', 'FontSize', 14)
text(1.05*x(2), 1.5*(y(2)+y(1))/2, '\Leftarrow Increasing Time', 'Rotation', 90, 'FontSize', 14)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 7 6])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the windowed correlation recovery
recoveryBlock = neuroCovWinBreaks(1:dayEnds(1), :);
recoveryBlock(recoveryBlock ==1) = NaN;
recovery = nanmean(recoveryBlock,1);
recoveryA = nanmean(windowTS(recovery, 51, 1, 'pad','boxcar')');
% recoveryB = nanmean(windowTS(recovery(2+1:end), 51, 1, 'pad','boxcar')');
% recovery = [recoveryA,recoveryB];
recovery = [recoveryA];
h(4) = figure(4);
cla;
plot(recovery)
hold on

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Windowed Correlation', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);


%%%%%%%%%%%%
%Save stuff to file
%Save images
% savefig(h, [sourceDir, figTitle, '.fig'])
% close(h)

%Save useful data
data = [];
data.annotFilelist = annotFilelist;
data.files = files;
data.dayEnds = dayEnds;
%data.lesionNum = lesionNum;
data.neuropowAbs = neuropowAbs;
data.neuropowAbs_winBreaks = neuropowAbs_winBreaks;
data.neuroCovWinBreaks = neuroCovWinBreaks;

save([sourceDir, figTitle, '.mat'], 'data')

clear all

