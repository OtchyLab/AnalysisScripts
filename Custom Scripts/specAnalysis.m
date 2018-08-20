%function specAnalysis
% 11/10/2015
% This script is for analying the gross structure of singing. It uses the koenigSpectral.m script to extract SAP features
% from intervals of sound defined by the annotation. There are a variety of setable options (at the beginning of code), the only input required is the
% name of the annotation file.

clear all

%Set data sources and analysis options
annotName = 'Grn229_150823_annotation_ALL.mat';
DOB = '07/04/2015';

startTime = '0:00';
endTime = '23:59';

%Generate all other file locations/names based on annotation name
mother = 'V:\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];
fileSource = [mother, char(sp{1}) filesep '20' char(sp{2}(1:2)) '-' char(sp{2}(3:4)) '-' char(sp{2}(5:6)) filesep];

outputName = [char(sp(1)) '_' char(sp(2))];
outputDir = ['C:\Users\Tim\Desktop\specAnalysis Output\' char(sp{1}) '\'];
if ~exist(outputDir)
    mkdir(outputDir) %create directory if it doesn't exist
end

%Calculate the relevants times
dataDay = sp{2};
dataDatenum = datenum(dataDay, 'yymmdd');
startDatenum = datenum([dataDay '_' startTime], 'yymmdd_HH:MM');
endDatenum = datenum([dataDay '_' endTime], 'yymmdd_HH:MM');
dobDatenum = datenum(DOB, 'mm/dd/yyyy');
dph = dataDatenum - dobDatenum;

%Load the annotation from file
load(annotLoc)

%Define constants for feature extraction
fs = 44150;
gain = 1;
hp = 500;
lp = 6500;

%Calculate constants for bandpass filtering
HP_fNorm = hp/(fs/2);
LP_fNorm = lp/(fs/2); 
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

%Establish output vars main vars 
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
    ins = repmat(keys(i),[numSyls,1]); %Speed up?
    wavName = [wavName; ins];
    
    %Unpack the acoustic features...
    bDoAcoustic = true;
    if (bDoAcoustic)
        %Calc the SAP features
        if strcmp(keys{i}(end), 'v')
            %It's a WAV
            features = koenigSpectral(filtfilt(BP_b,BP_a,audioread([fileSource keys{i}]).*gain), fs);
        else
            %It's a DAT
            [chans, ~] = getChannels([fileSource keys{i}]);
            features = koenigSpectral(filtfilt(BP_b,BP_a,chans(1,:).*gain), fs);
        end
        
        %Parse by syllables
        ins = zeros(1,numSyls);
        for j = 1:numSyls
            s = max([1, elements{i}.segFileStartTimes(j)*1000]);
            e = min([length(features.Entropy) elements{i}.segFileEndTimes(j)*1000]);
            snip = floor(s):ceil(e);
            ins(j) = mean(features.Entropy(snip)); %Change here to modify the features extracted
        end
        sylFeat = [sylFeat; ins'];
    end
end
sylDur = sylDur*1000; %Convert to ms

%Create selection indices
sylInx = (sylType>=1 & sylType<101) | sylType==103; %Include anything but call and 'unlabeled'
timeInx = sylAbsStartTime>=startDatenum & sylAbsStartTime<=endDatenum;

% %Calculate duration distributions for the three classes
minval = 0;
maxval = 350;
binSize = 1.25; %bin size in ms
[durBins,durCnts] = epdf_cbins(sylDur(sylInx & timeInx),binSize,minval,maxval);

%Plot the duration distributions on the same axis
h(1) = figure(1); clf

subplot(2,1,1)
plot(durBins,durCnts, 'LineWidth', 2)
%Format the figure
axis tight
xlim([0 350])
ylim([0,max(durCnts)*1.1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 100:100:350, 'LineWidth', 3, 'FontSize', 16)
xlabel('Syl Dur (ms)', 'FontSize', 16)
ylabel('P(t)', 'FontSize', 16)
%set(gcf, 'Units', 'Inches', 'Position', [0 0 5 5])
title([sp{1} ' ' sp{2} ' ' num2str(dph) 'dph'], 'FontSize', 16);

%Plot heatmaps for the durations and feature stats
if bDoAcoustic
    subplot(2,1,2)
    N = ndhist(sylDur(sylInx & timeInx), sylFeat(sylInx & timeInx), 'edgesx', 0:1.25:350, 'edgesy', -4:0.025:0, 'prob', 'axis', [-4 0 0 350]);
    sc = caxis; colormap(jet)
    %Format the figure
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 100:100:350, 'YTick', -4:2:0, 'LineWidth', 3, 'FontSize', 16)
    xlabel('Syl Dur (ms)', 'FontSize', 16)
    ylabel('Log(Entr)', 'FontSize', 16)
    set(gcf, 'Units', 'Inches', 'Position', [0 0 5 6])
end

%Save all data
save([outputDir, outputName, '_fullDaySpecAnalysis.mat']);

%Save figures
savefig(h, [outputDir, outputName, '_fullDaySpecAnalysis.fig']);
