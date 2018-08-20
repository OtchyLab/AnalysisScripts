function makeWaterfallHiRes
% The script creates a waterfall plot of MUA recordings that have been aligned with StretchEM.
% It solely looks at the neurtal recording data... use LongSong to look at the song data for this experiment.
% 
% This version removes the unnecessary plots (euclidean distance, etc) from the original, updates the algorithm for finding
% critical points (from indicies to time/date based), and adds in some code for removing noisy recordings.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Lesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn162\';
% annotFilelist = {'Grn162_150126_dataset ch4.mat';...
%                         'Grn162_150127_dataset ch4.mat'};
% 
% figTitle = 'Grn162 HVC Activity (Ch4)';
% lesionNum = 739;%

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn156\';
% annotFilelist = {'Grn156_150213_dataset ch2.mat';...
%                         'Grn156_150216_dataset ch2.mat'};
% 
% figTitle = 'Grn156 HVC Activity (Ch2)';
% lesionNum = 187; %

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn145\';
% annotFilelist = {'Grn145_150220_dataset ch4.mat';...
%                         'Grn145_150221_dataset ch4.mat';...
%                         'Grn145_150222_dataset ch4.mat';...
%                         'Grn145_150224_dataset ch4.mat'};
% 
% figTitle = 'Grn145 HVC Activity (Ch4)';
% lesionNum = 881; %881 as of 2/24/2015
 
% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn043\';
% annotFilelist = {'Grn043_140508_dataset ch4.mat';...
%                         'Grn043_140511_dataset ch4.mat'};
% 
% figTitle = 'Grn043 HVC Activity (Ch4)';
% lesionNum = 158; %158 as of 1/19/2015

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn051\';
% annotFilelist = {'Grn051_140528_dataset ch2.mat';...
%                         'Grn051_140531_dataset ch2.mat'};
% 
% figTitle = 'Grn051 HVC Activity (Ch2)';
% lesionNum = 117; %117 as of 1/20/2015

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn052\';
% annotFilelist = {'Grn052_140614_dataset ch4.mat';...
%                         'Grn052_140617_dataset ch4.mat'};
% 
% figTitle = 'Grn052 HVC Activity (Ch4)';
% lesionNum = 200; %680 as of 1/6/2015

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn094\';
% annotFilelist = {'Grn094_140619_dataset ch4.mat';...
%                         'Grn094_140623_dataset ch4.mat'};
% 
% figTitle = 'Grn094 HVC Activity (Ch4)';
% lesionNum = 845; %845 as of 1/20/2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn046\';
% annotFilelist = {'Grn046_140509_dataset ch3.mat';...
%                         'Grn046_140510_dataset ch3.mat';...
%                         'Grn046_140513_dataset ch3.mat';...
%                         'Grn046_140514_dataset ch3.mat'};
% 
% figTitle = 'Control046 HVC Activity (Ch3)';
% critTime = datenum(2014, 5, 09, 12, 30, 00); %5/15/2014 @ 11:32

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn121\';
% annotFilelist = {'Grn121_140831_dataset ch3.mat';...
%                         'Grn121_140902_dataset ch3.mat';...
%                         'Grn121_140904_dataset ch3.mat';...
%                         'Grn121_140905_dataset ch3.mat'};
% 
% figTitle = 'Control121 HVC Activity (Ch3)';
% critTime = datenum(2014, 8, 31, 20, 00, 00); %9/6/2014 2:05pm

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn141\';
% annotFilelist = {'Grn141_141207_dataset ch3.mat';...
%                         'Grn141_141208_dataset ch3.mat';...
%                         'Grn141_141209_dataset ch3.mat';...
%                         'Grn141_141210_dataset ch3.mat'};
% 
% figTitle = 'Control141 HVC Activity (Ch3)';
% critTime = datenum(2014, 12, 07, 14, 00, 00); %12/10/14 11:55am

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Complete Lesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn046\';
% annotFilelist = {'Grn046_140515_dataset ch3.mat';...
%                         'Grn046_140516_dataset ch3.mat';...
%                         'Grn046_140517_dataset ch3.mat';...
%                         'Grn046_140518_dataset ch3.mat'};
% 
% figTitle = 'Grn046 HVC Activity (Ch3)';
% % lesionNum = 385; %385 as of 1/6/2015
% critTime = datenum(2014, 5, 15, 11, 32, 00); %5/15/2014 @ 11:32

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn121\';
% annotFilelist = {'Grn121_140906_dataset ch3.mat';...
%                         'Grn121_140907_dataset ch3.mat';...
%                         'Grn121_140908_dataset ch3.mat';...
%                         'Grn121_140909_dataset ch3.mat'};
% 
% figTitle = 'Grn121 HVC Activity (Ch3)';
% % lesionNum = 1411; %1410 as of 1/13/2015
% critTime = datenum(2014, 9, 6, 14, 5, 00); %9/6/2014 2:05pm

sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn141\';
annotFilelist = {'Grn141_141210_dataset ch4.mat';...
                        'Grn141_141211_dataset ch4.mat';...
                        'Grn141_141212_dataset ch4.mat';...
                        'Grn141_141214_dataset ch4.mat'};

figTitle = 'Grn141 HVC Activity (Ch4)';
% lesionNum = 680; %680 as of 1/6/2015
critTime = datenum(2014, 12, 10, 11, 55, 00); %12/10/14 11:55am

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn186\';
% annotFilelist = {'Grn186_150424_dataset ch4.mat';...
%                        % 'Grn186_150425_dataset ch4.mat';...
%                         'Grn186_150426_dataset ch4.mat';...
%                         %'Grn186_150427_dataset ch4.mat';...
%                         'Grn186_150428_dataset ch4.mat';...
%                         'Grn186_150429_dataset ch4.mat'};
% figTitle = 'Grn186 HVC Activity (Ch4)';
% critTime = datenum(2015, 4, 24, 11, 50, 00); %

% Predefine variables
files = []; neuropow = []; neuropowAbs = []; dayEnds = []; dates = [];
                    
% Cycle through the files, loading into memory what is needed and deleting the rest
numAnnots = length(annotFilelist);
trimLength = 10; %in ms
trim = floor((trimLength/1000)*44150);
for i = 1:numAnnots
    % Clean load from drive
    data = []; filenames = [];
    load([sourceDir, annotFilelist{i}], 'data', 'filenames');
    
    %Filter out noisy recordings
    meanrawPow = mean(data.neuroPower_abs, 2);

    %shortwave
    if length(meanrawPow) > 100
        mean200 = nanmean(windowTS(meanrawPow, 201,1,'pad','boxcar')');
        std200 = nanstd(windowTS(meanrawPow, 201,1,'pad','boxcar')');
        thresh = mean200 + .75*std200;
    else
        mean200 = nanmean(windowTS(meanrawPow, 51,1,'pad','boxcar')');
        std200 = nanstd(windowTS(meanrawPow, 51,1,'pad','boxcar')');
        thresh = mean200 + .75*std200;
    end
    
    cleanIdx = meanrawPow<thresh';
    
    %Add relevant data to the running stacks
    files = [files, filenames(cleanIdx)];
    neuropowAbs = [neuropowAbs; data.neuroPower_abs(cleanIdx, trim:(end-trim))];
    
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
tmp = ind(timeVect<critTime);
lesionNum = tmp(end);
dayEnds = sort([dayEnds, lesionNum]);

%Calculate the mean power
meanPow = mean(neuropowAbs, 2);

% Generate covariance matrices using windowing
winSize = 25; winStep = 1;
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
imagesc(neuropowAbs_winBreaks, [0,.2]); % waterfall
hold on

%Format the subplot
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 18)
ylabel('Rendition', 'FontSize', 18)
title(figTitle, 'FontSize', 18);

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
set(gcf, 'Units', 'Inches', 'Position', [0 0 5 6]);


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
recoveryBlock = neuroCovWinBreaks(1:lesionNum, :);
recoveryBlock(recoveryBlock ==1) = NaN;
recovery = nanmean(recoveryBlock,1);
recoveryA = nanmean(windowTS(recovery(1:lesionNum), 51, 1, 'pad','boxcar')');
recoveryB = nanmean(windowTS(recovery(lesionNum+1:end), 51, 1, 'pad','boxcar')');
recovery = [recoveryA,recoveryB];
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
savefig(h, [sourceDir, figTitle, '.fig'])
% close(h)

%Save useful data
data = [];
data.annotFilelist = annotFilelist;
data.files = files;
data.dayEnds = dayEnds;
data.lesionNum = lesionNum;
data.neuropowAbs = neuropowAbs;
data.neuropowAbs_winBreaks = neuropowAbs_winBreaks;
data.neuroCovWinBreaks = neuroCovWinBreaks;

save([sourceDir, figTitle, '.mat'], 'data')

clear all

