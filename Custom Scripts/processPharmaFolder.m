function processPharmaFolder

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\Pharma\');

%Get all data files that match the criterion
crit = '*pharmaDataset.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially process the datasets
mcSet = [];
for i = 1:length(filelist)
    dataName = filelist(i).name;
    
    %Generate summary stats
%          pharmaSumStats(dataName, dataFolder);
    
    %Fix the datasets
    %     fixNumbering(dataName, dataFolder);
    %    update4Counts(dataName, dataFolder);
    %     includeMC(dataName, dataFolder);
    includeMC2(dataName, dataFolder);
    
    %     mcSet = getMC(dataName, dataFolder, mcSet);
    
end
% folder = 'C:\Users\Tim\Desktop\Nif Project Figures\Motif Continuity Analysis';
% save([folder, filesep, 'Control089 MCanal.mat'], 'mcSet')

display(['Finished with folder: ' dataFolder])


function mcSet = getMC(dataName, dataFolder,mcSet)
%Load data from file
load([dataFolder filesep dataName])
close all

mcSet = [mcSet;offset, meanMC];

function includeMC2(dataName, dataFolder2)
%Adds the motif continuity statistic to the pharamcologically lesioned datasets. Code is taken (and adapted) from the code
%for electrolytically lesioned data

%Load data from file
load([dataFolder2 filesep dataName])

% if exist('MC') ~=0 %check if present
%%%%%%%%%%%%%%%%%%%%%
%Set up sequencing and motif bands
%%%%%%%%%%%%%%%%%%%%%
% motif = [1,2,3,4,5,6];%Grn010
% motif = [1,2,3,4];%Grn011
motif = [3,4,2];%Grn089
% motif = [1,2,3,4];%Grn091
% motif = [2,3,4];%Pur935
date = char(sp{2});

mornEnd = '14:00:00';
nightStart = '19:00:00';

%Number of motifs within a sequencing window
winSize = 201;

%Reset Calculated Variables
sylTrans = [];
sylTransTime = [];
sylSeq = [];
sylSeqTime = [];
bands = [];
bandTime = [];
transMats = [];

%It's some other day; keep elements intact
b = [1, numRends];

%Setup the windowing for analysis
for j = 1:size(b,1)
    %Calculate sequencing primitives
    [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));
    
    %Push to permenant arrays
    sylSeq = [sylSeq; segTypeExt];
    sylSeqTime = [sylSeqTime; timestamp];
    sylTrans = [sylTrans; allSeq];
    sylTransTime = [sylTransTime; allSeqTime'];
end

%For summary stats
mornSeqInx = sylTransTime < datenum([date mornEnd],'yymmddHH:MM:SS');
nightSeqInx = sylTransTime > datenum([date nightStart],'yymmddHH:MM:SS');

%Motif Continuity -- probability that once you start the sequence, you'll finish it
mornMC = []; nightMC = [];
%Morning
tempArray = sylTrans(mornSeqInx,1)';
for j = 1:length(motif)
    if j == 1
        %This condition deals with variable repeats/stutters of the opening syllable
        pntrs = strfind(tempArray, motif(1));
        if ~isempty(pntrs) && pntrs(end) == length(tempArray)
            pntrs = pntrs(1:end-1);
        end
        pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
        mornMC(j) = length(pntrs_refined);
    else
        mornMC(j) = length(strfind(tempArray, motif(1:j)));
    end
end
jammed = length(strfind(tempArray, [motif(1:(end-1)), 101]));
if jammed~=0
    display('jammed syllables found')
    mornMC = [mornMC(1:end-1)-jammed, mornMC(end)];
end
mornMC = mornMC./mornMC(1);

%Night
tempArray = sylTrans(nightSeqInx,1)';
for j = 1:length(motif)
    if j == 1
        %This condition deals with variable repeats/stutters of the opening syllable
        pntrs = strfind(tempArray, motif(1));
        if ~isempty(pntrs) && pntrs(end) == length(tempArray)
            pntrs = pntrs(1:end-1);
        end
        pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
        nightMC(j) = length(pntrs_refined);
    else
        nightMC(j) = length(strfind(tempArray, motif(1:j)));
    end
end
jammed = length(strfind(tempArray, [motif(1:(end-1)), 101]));
if jammed~=0
    display('jammed syllables found')
    nightMC = [nightMC(1:end-1)-jammed, nightMC(end)];
end
nightMC = nightMC./nightMC(1);


meanMC = [mornMC(end), nightMC(end)];

%Save all data
save([dataFolder2 filesep dataName]);

%Save figures
%     savefig(h, [outputDir, outputName, '.fig']);
% end

%Close any open windows
close all



function includeMC(dataName, dataFolder2)
%Adds the motif continuity statistic to the pharamcologically lesioned datasets. Code is taken (and adapted) from the code
%for electrolytically lesioned data

%Load data from file
load([dataFolder2 filesep dataName])

% if exist('MC') ~=0 %check if present
%%%%%%%%%%%%%%%%%%%%%
%Set up sequencing and motif bands
%%%%%%%%%%%%%%%%%%%%%
% motif = [1,2,3,4,5,6];%Grn010
% motif = [1,2,3,4];%Grn011
motif = [3,4,2];%Grn089
% motif = [1,2,3,4];%Grn091
% motif = [2,3,4];%Pur935
date = char(sp{2});

mornEnd = '14:00:00';
nightStart = '19:00:00';

%Number of motifs within a sequencing window
winSize = 201;

%Reset Calculated Variables
sylTrans = [];
sylTransTime = [];
sylSeq = [];
sylSeqTime = [];
bands = [];
bandTime = [];
transMats = [];

%It's some other day; keep elements intact
b = [1, numRends];

%Setup the windowing for analysis
for j = 1:size(b,1)
    %Calculate sequencing primitives
    [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));
    
    %Push to permenant arrays
    sylSeq = [sylSeq; segTypeExt];
    sylSeqTime = [sylSeqTime; timestamp];
    sylTrans = [sylTrans; allSeq];
    sylTransTime = [sylTransTime; allSeqTime'];
end

%For summary stats
mornSeqInx = sylTransTime < datenum([date mornEnd],'yymmddHH:MM:SS');
nightSeqInx = sylTransTime > datenum([date nightStart],'yymmddHH:MM:SS');

%Motif Continuity -- probability that once you start the sequence, you'll finish it
mornMC = []; nightMC = [];
%Morning
tempArray = sylTrans(mornSeqInx,1)';
for j = 1:length(motif)
    if j == 1
        %This condition deals with variable repeats/stutters of the opening syllable
        pntrs = strfind(tempArray, motif(1));
        if ~isempty(pntrs) && pntrs(end) == length(tempArray)
            pntrs = pntrs(1:end-1);
        end
        pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
        mornMC(j) = length(pntrs_refined);
    else
        mornMC(j) = length(strfind(tempArray, motif(1:j)));
    end
end
jammed = length(strfind(tempArray, [motif(1:(end-1)), 101]));
if jammed~=0
    display('jammed syllables found')
    mornMC = [mornMC(1:end-1)-jammed, mornMC(end)];
end
mornMC = mornMC./mornMC(1);

%Night
tempArray = sylTrans(nightSeqInx,1)';
for j = 1:length(motif)
    if j == 1
        %This condition deals with variable repeats/stutters of the opening syllable
        pntrs = strfind(tempArray, motif(1));
        if ~isempty(pntrs) && pntrs(end) == length(tempArray)
            pntrs = pntrs(1:end-1);
        end
        pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
        nightMC(j) = length(pntrs_refined);
    else
        nightMC(j) = length(strfind(tempArray, motif(1:j)));
    end
end
jammed = length(strfind(tempArray, [motif(1:(end-1)), 101]));
if jammed~=0
    display('jammed syllables found')
    nightMC = [nightMC(1:end-1)-jammed, nightMC(end)];
end
nightMC = nightMC./nightMC(1);


meanMC = [mornMC(end), nightMC(end)];

%Save all data
save([dataFolder2 filesep dataName]);

%Save figures
%     savefig(h, [outputDir, outputName, '.fig']);
% end

%Close any open windows
close all

function update4Counts(dataName, dataFolder)

%Load from file
load([dataFolder filesep dataName])

%Do what needs to be updated
%%%%%%%%%%%%%%%\
minval = 0;
maxval = 350;
binSize = 1.25; %bin size in ms
[pdfX,pdfY] = epdf_cbins(sylDur(sylInx),binSize,minval,maxval);

%Plot the duration distributions on the same axis
h(1) = figure(1); clf; hold on
plot(pdfX,pdfY, 'LineWidth', 2)

%Format the figure
axis tight
xlim([0 350])
ylim([0,max(pdfY)*1.1])
box off
set(gca, 'TickDir', 'out', 'LineWidth', 3, 'FontSize', 16)
xlabel('Syl Dur (ms)', 'FontSize', 16)
ylabel('P(t)', 'FontSize', 16)

set(gcf, 'Units', 'Inches', 'Position', [0 0 4 3])

%Plot scatter for the durations
h(2) = figure(2); clf
scatter((sylAbsStartTime(sylInx)-dayDatenum)*24, sylDur(sylInx), '.')
hold on
axis tight
ylim([0, 350])
box off
set(gca, 'TickDir', 'out', 'LineWidth', 3, 'FontSize', 16)
xlabel('Time of Day (hrs)', 'FontSize', 16)
ylabel('Syl Dur (ms)', 'FontSize', 16)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 6 3])

%Plot heatmaps for the durations and feature stats
if bDoAcoustic
    h(3) = figure(3); clf
    
    N = ndhist(sylFeat(sylInx), sylDur(sylInx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'axis', [-4 0 0 350]);
    colormap(jet)
    box off
    set(gca, 'TickDir', 'out', 'YTick', 100:100:350, 'LineWidth', 3, 'FontSize', 16)
    ylabel('Syl Dur (ms)', 'FontSize', 16)
    xlabel('Mean log(Entropy)', 'FontSize', 16)
    title([sp{1} ' ' sp{2} ' ' conditions{conIdx}], 'FontSize', 16);

    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 4 4])
end

%Plot heatmap of durations
% Filter and mapping parameters
h(4) = figure(4); clf
[matrix, upperScale] = makeDurHeatmap(sylDur(sylInx));
imagesc(matrix, [0 upperScale]);
hold on
y = ylim;
colormap(jet)
box off
set(gca, 'TickDir', 'out', 'YTick', 1000:1000:size(matrix,1), 'LineWidth', 3, 'FontSize', 16)
xlabel('Syllable Durations (ms)', 'FontSize', 16)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 6])

%Save the output
save([dataFolder filesep dataName])

close all

function fixNumbering(dataName, dataFolder)
%Corrects the dataset indexing now that I've added Lesion-4/5d to the dataset.

%Load from file
load([dataFolder filesep dataName])
close all

if length(conditions) < 11 %This was the major change
    conditions = [{'Lesion - 5d'}, {'Lesion - 4d'}, {'Lesion - 3d'},{'Lesion - 2d'},{'Lesion - 1d'}, {'Lesion Day'}, {'Lesion + 1d'},{'Lesion + 2d'},{'Lesion + 3d'},{'Lesion + 4d'},{'Lesion + 5d'}];
    conIdx = offset+6;
    
    save([dataFolder filesep dataName])
end