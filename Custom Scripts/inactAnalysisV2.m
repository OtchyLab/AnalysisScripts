%function inactAnalysisV2
%3/2/2015
%This script is for analying changes in song before, during, and after inactivations/injections. It corrects some
%of the binning and plotting issues found in the original script (i.e., inactAnalysis.m).  

clear all

%Set data sources and injection time
annotName = 'Grn117_140820_annotation.mat';
inactStartTime = '12:45';
expType = 'inact';

%Define main vars
inactDur = 1; %time post injection in hours
postWash = 6; %time post injection in hours

%Generate all other file locations/names based on annotation name
mother = 'V:\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];
fileSource = [mother, char(sp{1}) filesep '20' char(sp{2}(1:2)) '-' char(sp{2}(3:4)) '-' char(sp{2}(5:6)) filesep];

outputDir = ['C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\' char(sp{1}) '\'];
outputName = [char(sp(1)) '_' char(sp(2))];

%Load the annotation from file
load(annotLoc)

%Define main vars for feature extraction
fs = 44150;
sylAbsStartTime = [];
sylDur = [];
sylType = [];
wavName = [];
sylFeat = [];
numRends = length(keys);
for i = 1:numRends
    %Unpack the time, sylType, and duration
    sylAbsStartTime = [sylAbsStartTime; elements{i}.segAbsStartTimes'];
    sylDur = [sylDur; [elements{i}.segFileEndTimes - elements{i}.segFileStartTimes]'];
    sylType = [sylType; elements{i}.segType];

    %Capture the filename of all syllables
    numSyls = length(elements{i}.segType);
    ins = cell(1,numSyls);
    for j = 1:numSyls
        ins{j} = keys{i};
    end
    wavName = [wavName; ins'];
    
    %Unpack the acoustic features...
    bDoAcoustic = true;
    if (bDoAcoustic)
        %Calc the SAP features
        if strcmp(keys{i}(end), 'v')
            %It's a WAV
            features = koenigSpectral(audioread([fileSource keys{i}]), fs);
        else
            %It's a DAT
            [chans, ~] = getChannels([fileSource keys{i}]);
            features = koenigSpectral(chans(1,:), fs);
        end

        %Parse by syllables
        ins = zeros(1,numSyls);
        for j = 1:numSyls
            s = max([1, elements{i}.segFileStartTimes(j)*1000]);
            e = min([length(features.Entropy) elements{i}.segFileEndTimes(j)*1000]);
            snip = floor(s):ceil(e);
            ins(j) = mean(features.Entropy(snip));
        end
        sylFeat = [sylFeat; ins'];
    end
    
end
sylDur = sylDur*1000; %Convert to ms

sp = regexp(annotName, '_', 'split');
inactStartDay = sp{2};
inactDatenum = datenum([inactStartDay '_' inactStartTime], 'yymmdd_HH:MM');

%Create indices
sylInx = (sylType>=1 & sylType<101) | sylType==103;
preInx = sylAbsStartTime<inactDatenum;
injInx = sylAbsStartTime>=inactDatenum & sylAbsStartTime<=(inactDatenum+(inactDur/24));
postInx = sylAbsStartTime>(inactDatenum+((postWash)/24));

%Calculate duration distributions for the three classes
minval = 0;
maxval = 350;
binSize = 1.25; %bin size in ms
[preX,preY] = epdf_cbins(sylDur(preInx & sylInx),binSize,minval,maxval);
[postX,postY] = epdf_cbins(sylDur(postInx & sylInx),binSize,minval,maxval);
[injX,injY] = epdf_cbins(sylDur(injInx & sylInx),binSize,minval,maxval);

%Plot the duration distributions on the same axis
h(1) = figure(1);
cla; hold on
plot(preX,preY, 'LineWidth', 2)
plot(postX,postY, 'LineWidth', 2)
plot(injX,injY, 'LineWidth', 2)

%Format the figure
axis tight
xlim([0 350])
ylim([0,max([preY, postY, injY])*1.1])
box off
set(gca, 'TickDir', 'out')
xlabel('Syl Dur (ms)', 'FontSize', 16)
ylabel('P(t)', 'FontSize', 16)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 5])
set(gca, 'LineWidth', 3, 'FontSize', 16)
title([sp{1} ' ' sp{2}], 'FontSize', 16);
legend({'pre'; 'post'; expType}); legend('boxoff')

%Plot scatter for the durations
h(2) = figure(2);
cla
scatter((sylAbsStartTime-inactDatenum)*24, sylDur, '.')
hold on
line([0, 0], [0, 350], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r')
xsl = xlim;
patch([0, inactDur, inactDur, 0], [0, 0, 350, 350], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([postWash, xsl(2), xsl(2), postWash], [0, 0, 350, 350], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([xsl(1), 0, 0, xsl(1)], [0, 0, 350, 350], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
axis tight
ylim([0, 350])
box off
set(gca, 'TickDir', 'out')
xlabel('Time (hrs)', 'FontSize', 16)
ylabel('Syl Dur (ms)', 'FontSize', 16)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 7 4])
set(gca, 'LineWidth', 3, 'FontSize', 16)
title([sp{1} ' ' sp{2}], 'FontSize', 16);


%Plot heatmaps for the durations and feature stats
if bDoAcoustic
    h(3) = figure(3); clf

    subplot(1,3,1)
    preN = ndhist(sylFeat(preInx & sylInx), sylDur(preInx & sylInx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'filter', 'axis', [-4 0 0 350]);
    sc = caxis; colormap(jet)
    box off
    set(gca, 'TickDir', 'out', 'YTick', 100:100:350)
    ylabel('Syl Dur (ms)', 'FontSize', 16)
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    legend('pre'); legend('boxoff')

    subplot(1,3,2)
    injN = ndhist(sylFeat(injInx & sylInx), sylDur(injInx & sylInx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'filter', 'axis', [-4 0 0 350]);
    caxis(sc); colormap(jet)
    box off
    set(gca, 'TickDir', 'out', 'YTickLabels', [], 'YTick', 100:100:350)
    xlabel('Mean log(Entropy)', 'FontSize', 16)
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    title([sp{1} ' ' sp{2}], 'FontSize', 16);
    legend(expType); legend('boxoff')

    subplot(1,3,3)
    postN = ndhist(sylFeat(postInx & sylInx), sylDur(postInx & sylInx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'filter', 'axis', [-4 0 0 350]);
    caxis(sc); colormap(jet)
    box off
    set(gca, 'TickDir', 'out', 'YTickLabels', [], 'YTick', 100:100:350)
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    legend('post'); legend('boxoff')

    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 8 5])
end

%Plot heatmap of durations
% Filter and mapping parameters
[matrix, upperScale] = makeDurHeatmap(sylDur(sylInx));
injMarkerPnt = find(sylAbsStartTime<=inactDatenum, 1, 'last');

h(4) = figure(4); clf
%The duration waterfall plot
subplot(1,5,1:4)
imagesc(matrix, [0 upperScale]);
hold on
line([0, 350], [injMarkerPnt, injMarkerPnt], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r')
y = ylim;
colormap(jet)
box off
set(gca, 'TickDir', 'out', 'YTick', 1000:1000:size(matrix,1))
xlabel('Syllable Durations (ms)', 'FontSize', 16)
title([sp{1} ' ' sp{2}], 'FontSize', 16);
set(gca, 'LineWidth', 3, 'FontSize', 16)

%The colorbar showing snip locations
prePatch =[find(preInx == true, 1, 'first'), find(preInx == true, 1, 'last')];
injPatch = [find(injInx == true, 1, 'first'), find(injInx == true, 1, 'last')];
postPatch = [find(postInx == true, 1, 'first'), find(postInx == true, 1, 'last')];

subplot(1,5,5); cla
hold on
patch([0,0,1,1], [prePatch fliplr(prePatch)], [0, 0.45, 0.74], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([0,0,1,1], [postPatch fliplr(postPatch)], [0.85, 0.33, 0.1], 'FaceAlpha', 1, 'EdgeColor', 'none')
patch([0,0,1,1], [injPatch fliplr(injPatch)], [0.93, 0.69, 0.13], 'FaceAlpha', 1, 'EdgeColor', 'none')
axis ij; axis tight
ylim(y); box off
set(gca, 'XTick', [], 'YTick', [])

set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7])

%Save all data
save([outputDir, outputName, '_inactDataset.mat']);

%Save figures
savefig(h, [outputDir, outputName, '.fig']);
