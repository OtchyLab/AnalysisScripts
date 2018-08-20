function lesionPartialEffect
%4-21-15 (TMO)
%The function of this script is to analyze the effect of ELECTROLYTIC lesions to Nif on HVC activity on a single day (as opposed to
%recovery over multiple days). It takes as input the aligned datasets from StretchEm.m
close all

%Per-bird data to update
annotName = 'Grn132_140918_dataset ch3.mat'; %23
lesionTime = datenum(2014, 9, 18, 13, 15, 00); %
lesionFrac = 0; %Percentage of Nif lesioned

%Setup for file locations
% mother = 'C:\Users\Tim\Desktop\Temp Nif Data\';
mother = 'E:\Active Datasets\Partial Nif Lesions\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];

%This is the amount to slice off of the fore- and aft sections of the waterfall to deal with edge effects of smoothing
trimLength = 10; %in ms
trim = floor((trimLength/1000)*44150);

% Clean load from drive
data = []; filenames = [];
load(annotLoc, 'data', 'filenames');

%Extra relevant data from the alignment file
files = filenames;
neuropowAbs = data.neuroPower_abs(:, trim:(end-trim));
clear('data', 'filenames')

%Extract the recording time from the filenames
fileTime = [];
for i = 1:length(files)
    fileTime(i) = getFileTime(files{i});
end

%Timing/sectioning control and counts
preTime = 2;    %hours before lesion to include
postTime = 1;   %hours of post-lesion singing to include
nightTimeS = 1;  %hours after lesion to begin
nightTimeE = 50;  %hours after lesion to end

lesionIdx = find(fileTime>lesionTime, '1', 'first'); %Time the bird starts singing after lesion
lesionIdx = lesionIdx(1);

%Clean the noise from the files
% [neuropowAbs, fileTime, noiseIdx] = idNoise(neuropowAbs,fileTime, lesionIdx);
[neuropowAbs, fileTime, noiseIdx] = idNoise(neuropowAbs,fileTime, lesionIdx);

%Recalc lesion index with clean files
lesionIdx = find(fileTime>lesionTime, '1', 'first'); %Time the bird starts singing after lesion
lesionIdx = lesionIdx(1);

%Make timing indices
postStartTime = fileTime(lesionIdx); 
preInx = fileTime>(lesionTime-(preTime/24)) & fileTime<=lesionTime;
postInx = fileTime>lesionTime & fileTime<(postStartTime+(postTime/24));
nightInx = fileTime>(lesionTime+(nightTimeS/24)) & fileTime<(lesionTime+(nightTimeE/24));

%Smooth the waterfall plot with a running window over renditions; don't cross the lesion boundary
winSize = 25; winStep = 1;                                                                   %windowing parameters
neuropowAbs_winBreaks = [];
sections = [1, lesionIdx-1; lesionIdx, size(neuropowAbs,1)];
for i = 1:size(sections,1)
    inxRend = sections(i,1):sections(i,2);                                                 %windowing index
    [indx_mat] = windowTS(inxRend,winSize,winStep,'pad','boxcar');
    temp = [];
    for j = 1:size(indx_mat,1)                                                                  %window it
        trimIndx = indx_mat(j,~isnan(indx_mat(j, :)));
        temp(j, :) = mean(neuropowAbs(trimIndx,:),1);
    end
    neuropowAbs_winBreaks = [neuropowAbs_winBreaks; temp];                                  %stack it
end

%Calculate the mean HVC power (per rendition)
meanPow = mean(neuropowAbs_winBreaks, 2);

% %Generate correlation matrix on the windowed data
% [neuroCovWin,~] = corrcoef(neuropowAbs_win');
% 
% %Calculate the recovery in correlation
% recovBlock = neuroCovWin(preInx, :);
% recovBlock(recovBlock ==1) = NaN;
% recovCorr = nanmean(recovBlock,1);
%
% %The colorbar showing snip locations
prePatch =[find(preInx == true, 1, 'first'), find(preInx == true, 1, 'last')];
postPatch = [find(postInx == true, 1, 'first'), find(postInx == true, 1, 'last')];
nightPatch = [find(nightInx == true, 1, 'first'), find(nightInx == true, 1, 'last')];

%New method for calculating the correlation recovery
%New Neural Dynamics recovery
%Assemble the matrix for correlations
origWin = 25;
winOffset = ceil(origWin/2);

pre = neuropowAbs_winBreaks(lesionIdx-winOffset, :);
post = neuropowAbs_winBreaks(lesionIdx+winOffset, :);
postEnd = neuropowAbs_winBreaks(end-winOffset, :);
neuralMat = [pre; post; postEnd];

%Take correlation of the pre-lesion with each subsequent rendition
neuroCorr = [];
for i = 2:size(neuralMat,1)
    neuroCorr(i) = corr(neuralMat(1,:)', neuralMat(i,:)');
end

%Caluclate the correlation between non-pverlapping pre-lesion windows
prepre = neuropowAbs_winBreaks(lesionIdx-2*origWin, :);
neuroCorr(1) = corr(neuralMat(1,:)', prepre');


%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting outputs
%%%%%%%%%%%%%%%%%%%%%%%%%
%Waterfall
h(1) = figure(3); clf;
subplot(1,5,1:4); cla;
imagesc(neuropowAbs_winBreaks); % waterfall
hold on
x = xlim; y = ylim;
    %Format the subplot
xticks = 100:100:(x(2)*1000/44150); xpos = round((xticks/1000)*44150);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', xpos, 'XTickLabels', xticks, 'LineWidth', 2, 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 18)
ylabel('Rendition', 'FontSize', 18)
%title(figTitle, 'FontSize', 18);

subplot(1,5,5); cla %This is the bar showing the breaks for each of the critical points
hold on
patch([0,0,1,1], [prePatch fliplr(prePatch)], [0, 0.45, 0.74], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([0,0,1,1], [postPatch fliplr(postPatch)], [0.85, 0.33, 0.1], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([0,0,1,1], [nightPatch fliplr(nightPatch)], [0.93, 0.69, 0.13], 'FaceAlpha', 1, 'EdgeColor', 'none')
axis ij; axis tight
ylim(y);
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7])

%Stats for the day
h(2) = figure(4); clf
subplot(2,1,1); cla
for i = 1:size(sections,1)
    plot(sections(i,1):sections(i,2), meanPow(sections(i,1):sections(i,2)),  '-k', 'LineWidth', 2.5); hold on
end
axis tight; x = xlim;
% xlim([0, 0.1])
ylabel('Mean HVC Activity (V^2)', 'FontSize', 8)
set(gca, 'Box', 'off',  'TickDir', 'out', 'YTick', [], 'LineWidth', 2, 'FontSize', 8)

% subplot(3,1,2); cla
% for i = 1:size(sections,1)
%     plot(sections(i,1):sections(i,2), recovCorr(sections(i,1):sections(i,2)),  '-k', 'LineWidth', 2.5); hold on
% end
% axis tight; xlim(x)
% % xlim([0, 0.1])
% ylabel('Corr to Pre-Lesion', 'FontSize', 8)
% set(gca, 'Box', 'off',  'TickDir', 'out', 'YTick', [], 'LineWidth', 2, 'FontSize', 8)

subplot(2,1,2); cla
hold on
patch([prePatch fliplr(prePatch)], [0,0,1,1], [0, 0.45, 0.74], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([postPatch fliplr(postPatch)], [0,0,1,1], [0.85, 0.33, 0.1], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([nightPatch fliplr(nightPatch)], [0,0,1,1], [0.93, 0.69, 0.13], 'FaceAlpha', 1, 'EdgeColor', 'none')
axis tight; xlim(x)
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [0 0 5 4])

%Plot the New correlation measure, unnromalized
h(3) = figure(5); clf
plot(1:length(neuroCorr), neuroCorr, '-ok')
xlabel('Time'); ylabel('Correlation');
xlim([0.5, 3.5])
set(gca, 'Box', 'off', 'XTick', 1:length(neuroCorr), 'XTickLabel', [{'PreLesion'}, {'PostLesion'}, {'Post Night'}])
set(gcf, 'Units', 'Inches', 'Position', [0 0 3 3])

%%%%%%%%%%%%%%%%%%%%%%
%Save workspace data
%%%%%%%%%%%%%%%%%%%%%%
outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\Lesion Size-Effect\';

%Save all data
save([outputLocation, sp{1} '-' sp{2}, '-ch', sp{3}(end-4), '_allData05192015.mat']);

%Save figures
% savefig(h, [outputLocation, sp{1} '-' sp{2}, '-ch', sp{3}(end-4), '_figures05192015.fig']);

%%%%%%%%%%%%%%%%%%%%%%
%Gather summary data for across-bird comparison
%%%%%%%%%%%%%%%%%%%%%%
output = [];
output.name = sp{1};
output.date = sp{2};
output.filename = annotName;
output.chan = str2num(sp{3}(end-4));
output.lesionFrac = lesionFrac;

    %HVC mean activity
    m = []; s = [];
m(1) = mean(meanPow(preInx));
m(2) = mean(meanPow(postInx));
m(3) = mean(meanPow(nightInx));
s(1) = std(meanPow(preInx));
s(2) = std(meanPow(postInx));
s(3) = std(meanPow(nightInx));
output.powMeans = m; output.powStd = s;

    %HVC recovery
%     m = []; s = [];
% m(1) = mean(recovCorr(preInx));
% m(2) = mean(recovCorr(postInx));
% m(3) = mean(recovCorr(nightInx));
% s(1) = std(recovCorr(preInx));
% s(2) = std(recovCorr(postInx));
% s(3) = std(recovCorr(nightInx));
% output.recoveryMeans = m; output.recoveryStd = s;

output.corrRecovery = neuroCorr;

%%%%%%%%%%%%%%%%%%%%%
% Save to output location
%%%%%%%%%%%%%%%%%%%%%

outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
outPerm = 'Partial Lesion Summary Data 05192015.mat';
m = exist([outputLocation, outPerm]);
if m ==2 
    %File already exists
    clear('partialSumStats');
    load([outputLocation, outPerm], 'partialSumStats')
    partialSumStats(end+1) = output;
else
    %No file yet created
    partialSumStats = output;
end

%Save the updated data to file
save([outputLocation, outPerm], 'partialSumStats')

%Write the noisy file names to a text file for removal later
crapFiles = files(noiseIdx);
fileID = fopen([outputLocation, sp{1} '-' sp{2}, '-ch', sp{3}(end-4), '_noisyFiles2.txt'],'w');
fprintf(fileID,'%1s \r\n','Name');
fprintf(fileID, '%1s \r\n', cell2mat(crapFiles)');
fclose(fileID);

display(['done with: ' annotName])

function [cleanWater, cleanFiles, noiseIdx] = idNoise(waterfall,filenames, lesionIdx)
% Cycle through the files, loading into memory what is needed and deleting the rest
neuropowAbs = [];
files = [];
trimLength = 10; %in ms
trim = floor((trimLength/1000)*44150);

%Filter out noisy recordings
meanrawPow = mean(waterfall, 2);

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
neuropowAbs = [neuropowAbs; waterfall(cleanIdx, trim:(end-trim))];

%Longwave
meantrimPow = mean(neuropowAbs,2);
mean1000 = nanmean(windowTS(meantrimPow, 501,1,'pad','boxcar')');
std1000 = nanstd(windowTS(meantrimPow,501,1,'pad','boxcar')');
thresh2 = mean1000 + 2*std1000;
cleanIdxLong = meantrimPow<thresh2';

%apply
cleanFiles = files(cleanIdxLong);
cleanWater = neuropowAbs(cleanIdxLong, :);
noiseIdx = ~ismember(cleanFiles,filenames);


function getBirdInfo(birdName)
%Complete lesion
annotName = 'Grn046_140515_dataset ch3.mat'; %234
lesionTime = datenum(2014, 5, 15, 11, 32, 00);
lesionFrac = 100; %Percentage of Nif lesioned

annotName = 'Grn121_140906_dataset ch3.mat'; %234
lesionTime = datenum(2014, 9, 6, 14, 5, 00);
lesionFrac = 100; %Percentage of Nif lesioned

annotName = 'Grn141_141210_dataset ch4.mat'; %234
lesionTime = datenum(2014, 12, 10, 11, 55, 00);
lesionFrac = 70; %Percentage of Nif lesioned

annotName = 'Grn186_150424_dataset ch4.mat'; %23
lesionTime = datenum(2015, 4, 24, 11, 50, 00);
lesionFrac = 70; %Percentage of Nif lesioned

%Partial lesion
annotName = 'Grn052_140614_dataset ch2.mat'; %234#
lesionTime = datenum(2014, 6, 14, 12, 56, 00);
lesionFrac = 10; %Percentage of Nif lesioned

annotName = 'Grn094_140619_dataset ch4.mat'; %234#
lesionTime = datenum(2014, 6, 19, 11, 44, 00);
lesionFrac = 10; %Percentage of Nif lesioned

annotName = 'Grn103_140721_dataset ch4.mat'; %234#
lesionTime = datenum(2014, 07, 21, 12, 48, 00);
lesionFrac = 30; %Percentage of Nif lesioned

annotName = 'Grn149_141226_dataset ch4.mat'; %234
lesionTime = datenum(2014, 12, 26, 13, 46, 00);
lesionFrac = 20; %Percentage of Nif lesioned

annotName = 'Grn156_150213_dataset ch4.mat'; %24
lesionTime = datenum(2015, 2, 13, 11, 12, 00); %
lesionFrac = 40; %Percentage of Nif lesioned

annotName = 'Grn165_150309_dataset ch4.mat'; %23
lesionTime = datenum(2015, 3, 9, 12, 21, 00); %
lesionFrac = 40; %Percentage of Nif lesioned



%No lesion
annotName = 'Grn051_140528_dataset ch4.mat'; %23
lesionTime = datenum(2014, 5, 28, 10, 30, 00); %
lesionFrac = 0; %Percentage of Nif lesioned

annotName = 'Grn132_140918_dataset ch4.mat'; %23
lesionTime = datenum(2014, 9, 18, 13, 15, 00); %
lesionFrac = 0; %Percentage of Nif lesioned

%Missed Lesion
annotName = 'Grn043_140508_dataset ch4.mat';
lesionTime = datenum(2014, 5, 8, 14, 12, 00); %
lesionFrac = 0; %Percentage of Nif lesioned

annotName = 'Grn132_140918_dataset ch4.mat';
lesionTime = datenum(2014, 9, 18, 13, 15, 00); %
lesionFrac = 0; %Percentage of Nif lesioned




