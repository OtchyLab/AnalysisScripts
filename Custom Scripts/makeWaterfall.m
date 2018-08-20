function [figTitle] = makeWaterfall(pntr, chan)
% The script creates a waterfall plot of MUA recordings that have been aligned with StretchEM

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn141\';
% annotFilelist = {'Grn141_141210_dataset ch2 300-6000Hz.mat';...
%                         'Grn141_141211_dataset ch2 300-6000Hz.mat';...
%                         'Grn141_141212_dataset ch2 300-6000Hz.mat';...
%                         'Grn141_141214_dataset ch2 300-6000Hz.mat'};
% 
% figTitle = 'Grn141 HVC Activity (Ch2)';
% lesionNum = 680; %680 as of 1/6/2015

sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn046\';
annotFilelist = {'Grn046_140515_dataset ch4.mat';...
                        'Grn046_140516_dataset ch4.mat';...
                        'Grn046_140517_dataset ch4.mat';...
                        'Grn046_140518_dataset ch4.mat'};

figTitle = 'Grn046 HVC Activity (Ch4)';
lesionNum = 385; %385 as of 1/6/2015

% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn121\';
% annotFilelist = {'Grn121_140906_dataset_ch4.mat';...
%                         'Grn121_140907_dataset_ch4.mat';...
%                         'Grn121_140908_dataset_ch4.mat';...
%                         'Grn121_140909_dataset_ch4.mat'};
% 
% figTitle = 'Grn121 HVC Activity (Ch4)';
% lesionNum = 1411; %1410 as of 1/13/2015
% 
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

% Set the source data
% sourceDir = 'C:\Users\Tim\Desktop\Temp Nif Data\Grn141\';
% lesionNum = 665;
% [annotFilelist, figTitle] = nameGen(pntr, chan);

% Predefine variables
files = []; neuropow = []; neuropowAbs = []; dayEnds = [];
                    
% Cycle through the files, loading into memory what is needed and deleting the rest
numAnnots = length(annotFilelist);
for i = 1:numAnnots
    % Clean load from drive
    data = []; filenames = [];
    load([sourceDir, annotFilelist{i}], 'data', 'filenames');
    
    %Add relevant data to the running stacks
    files = [files, filenames];
%    neuropow = [neuropow; data.neuroPower(:, 10:(end-10))];
    neuropowAbs = [neuropowAbs; data.neuroPower_abs(:, 10:(end-10))];
    
    dayEnds(i) = length(files);
end

% Add lesion marker
dayEnds = sort([dayEnds, lesionNum]);

%Calculate the mean power
meanPow = mean(neuropowAbs, 2);

%Generate covariance matrices
[neuroCov,neuroP] = corrcoef(neuropowAbs');

%Generate covariance matrices using windowing
winSize = 25; winStep = 1; inxRend = 1:size(neuropowAbs, 1);
[indx_mat] = windowTS(inxRend,winSize,winStep,'pad');
neuropowAbs_win = [];
for i = 1:size(indx_mat,1)
    trimIndx = indx_mat(i,~isnan(indx_mat(i, :)));
    neuropowAbs_win(i, :) = mean(neuropowAbs(trimIndx,:),1);
end

[neuroCovWin,neuroPWin] = corrcoef(neuropowAbs_win');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = [];

%Plot the waterfall image
h(1) = figure(1);
cla;
imagesc(neuropowAbs, [0,.3]); % waterfall
hold on

% Day markers
x = xlim;
t = [0, dayEnds];
for i = 2:length(t)
    plot([1.03*x(2), 1.03*x(2)], [(t(i-1)+0.5), t(i)-.5],'LineWidth', 4, 'clip', 'off')
end
% for i = 1:length(dayEnds)
%     plot(1.03*x(2), dayEnds(i), '<k', 'MarkerSize', 10, 'clip', 'off')
% end

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Time (ms)', 'FontSize', 18)
ylabel('Rendition', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 6])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

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

% Plot the mean power per day
h(3) = figure(3);
cla;
hold on

t = [0, dayEnds];
space = 0.25;
for i = 2:length(t)
    rendSnips = (t(i-1)+1):t(i);
    plot(mean(neuropowAbs(rendSnips,:), 1) - space*(i-2))
end
x = xlim;
y = ylim;

%Format the image
box off
set(gca, 'TickDir', 'out','YTickLabels', [])
xlabel('Time (ms)', 'FontSize', 18)
ylabel('Mean HVC Activity (V^2)', 'FontSize', 18)
text(1.05*x(2), 2.25*(y(2)+y(1))/2, '\Leftarrow Increasing Time', 'Rotation', 90, 'FontSize', 16)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the rendition by rendition correlation
h(4) = figure(4);
cla;

%Show data on axes
imagesc(neuroCov, [0, 1]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the windowed correlation
h(5) = figure(5);
cla;

%Show data on axes
imagesc(neuroCovWin, [0, 1]);
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


%%%%%%%%%%%%%%%%%%%%%
% Plot the correlation recovery
recoveryBlock = neuroCov(1:lesionNum, :);
recoveryBlock(recoveryBlock ==1) = NaN;
recovery = nanmean(recoveryBlock,1);
h(6) = figure(6);
cla;
plot(recovery)
hold on

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Correlation', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the windowed correlation recovery
recoveryBlock = neuroCovWin(1:lesionNum, :);
recoveryBlock(recoveryBlock ==1) = NaN;
recovery = nanmean(recoveryBlock,1);
h(7) = figure(7);
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
% Plot the euclidean distance matrix
euclidDist = squareform(pdist(neuropowAbs,'euclidean'));
h(8) = figure(8);
cla;
imagesc((euclidDist), [0, 1.5]);
hold on
for i = 2:length(t)
    plot([x(1) - .03*x(2), x(1) - .03*x(2)], [(t(i-1)+8), t(i)-8],'LineWidth', 4, 'clip', 'off')
end

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Euclidean Distance', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%
% Plot the windowed euclidean distance matrix
euclidDist_win = squareform(pdist(neuropowAbs_win,'euclidean'));
h(9) = figure(9);
cla;
imagesc((euclidDist_win), [0, 1.5]);
hold on
for i = 2:length(t)
    plot([x(1) - .03*x(2), x(1) - .03*x(2)], [(t(i-1)+8), t(i)-8],'LineWidth', 4, 'clip', 'off')
end

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Windowed Euclidean Distance', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%
% Plot the  euclid distance recovery
recoveryBlock = euclidDist(1:lesionNum, :);
recoveryBlock(recoveryBlock ==0) = NaN;
recovery = nanmean(recoveryBlock,1);
h(10) = figure(10);
cla;
plot(recovery)
hold on

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Euclidean Distance', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

%%%%%%%%%%%%%%%%%%%%%%%
% Plot the  euclid distance recovery
recoveryBlock = euclidDist_win(1:lesionNum, :);
recoveryBlock(recoveryBlock ==0) = NaN;
recovery = nanmean(recoveryBlock,1);
h(11) = figure(11);
cla;
plot(recovery)
hold on

%Format the image
box off
set(gca, 'TickDir', 'out')
xlabel('Renditions', 'FontSize', 18)
ylabel('Windowed Euclidean Distance', 'FontSize', 18)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 4])
set(gca, 'LineWidth', 2, 'FontSize', 14)
title(figTitle, 'FontSize', 18);

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
data.neuropowAbs_win = neuropowAbs_win;
data.neuroCov = neuroCov;
data.neuroCovWin = neuroCovWin;
data.euclidDist = euclidDist;
data.euclidDist_win = euclidDist_win;

save([sourceDir, figTitle, '.mat'], 'data')

clear all

