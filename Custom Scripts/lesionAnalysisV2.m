%3/28/2015
%This script is for analying 1D and 2D pdfs for song structure -- mainly for the PHARMACOLOGICAL
%Nif lesions and is generalized to collect summary info over many days.. It corrects some
%of the binning and plotting issues found in the original script (i.e., inactAnalysis.m)  
%
%Updated 5/11/2015 to double bin sizes on histograms


clear all; close all

%Set data sources and injection time
annotName = 'Pur935_131118_annotation.mat';
lesionDate = '2013-11-13';

%Generate all other file locations/names based on annotation name
mother = 'V:\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];
fileSource = [mother, char(sp{1}) filesep '20' char(sp{2}(1:2)) '-' char(sp{2}(3:4)) '-' char(sp{2}(5:6)) filesep];

outputDir = ['C:\Users\Tim\Desktop\Nif Project Figures\Pharma\' char(sp{1}) '\'];
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


%Determine which day (relative to lesion) this is
sp = regexp(annotName, '_', 'split');
Day = sp{2};
dayDatenum = datenum(Day, 'yymmdd');
lesionDatenum = datenum(lesionDate, 'yyyy-mm-dd');
offset = dayDatenum - lesionDatenum;
conditions = [{'Lesion - 5d'}, {'Lesion - 4d'}, {'Lesion - 3d'},{'Lesion - 2d'},{'Lesion - 1d'}, {'Lesion Day'}, {'Lesion + 1d'},{'Lesion + 2d'},{'Lesion + 3d'},{'Lesion + 4d'},{'Lesion + 5d'},{'Lesion + 6d'},{'Lesion + 7d'},{'Lesion + 8d'}];
conIdx = offset+6;

%Create syllable index
sylInx = (sylType>=1 & sylType<101) | sylType==103;

%Calculate duration distributions for the three classes
% minval = floor(min(sylDur));
% maxval = ceil(max(sylDur));
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
%title([sp{1} ' ' sp{2} ' ' conditions{conIdx}], 'FontSize', 16);

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 4 3])

%Plot scatter for the durations
h(2) = figure(2); clf
scatter((sylAbsStartTime(sylInx)-dayDatenum)*24, sylDur(sylInx), '.')
hold on
xsl = xlim;
axis tight
ylim([0, 350])
box off
set(gca, 'TickDir', 'out', 'LineWidth', 3, 'FontSize', 16)
xlabel('Time of Day (hrs)', 'FontSize', 16)
ylabel('Syl Dur (ms)', 'FontSize', 16)
% title([sp{1} ' ' sp{2} ' ' conditions{conIdx}], 'FontSize', 16);

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
%title([sp{1} ' ' sp{2} ' ' conditions{conIdx}], 'FontSize', 16);

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 6])

% %%%%%%%%%%%%%%%%%%%%%
% %Set up sequencing and motif bands
% %%%%%%%%%%%%%%%%%%%%%
% motif = [1,2,3,4,5,6];
% date = char(sp{2});
% 
% mornEnd = '12:00:00';
% nightStart = '19:00:00';
% 
% %Number of motifs within a sequencing window
% winSize = 201;
% 
% %Reset Calculated Variables
% sylTrans = [];
% sylTransTime = [];
% sylSeq = [];
% sylSeqTime = [];
% bands = [];
% bandTime = [];
% transMats = [];
% 
% %It's some other day; keep elements intact
% b = [1, numRends];
% 
% %Setup the windowing for analysis
% for j = 1:size(b,1)
%     %Calculate sequencing primitives
%     [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));
% 
%     %Find the starts to all motifs
%     Indx = strfind(allSeq(:,1)', motif);
% 
%     %Generate sequence start and ends
%     try
%        indWin = windowTS(Indx, winSize, 1, 'pad', 'boxcar');
%     catch
%         indWin = Indx;
%     end
%     motifBands = [indWin(:,1), indWin(:,end)+length(motif)-1];
%     motifBands(isnan(motifBands(:,1)),1) = Indx(1);
%     motifBands(isnan(motifBands(:,2)),2) = Indx(end)+length(motif)-1;
%     motifBands(motifBands<1) = 1;
%     motifBands(motifBands>size(allSeq,1)) = size(allSeq,1);
%     bands = [bands;motifBands+size(sylTrans,1)];
% 
%     %Locate the center time for each band
%     try
%         t = allSeqTime(indWin(:,ceil(winSize/2)));
%     catch
%         t = allSeqTime(indWin(ceil(length(indWin)/2)));
%     end
%     bandTime = [bandTime; t'];
% 
%     %Push to permenant arrays
%     sylSeq = [sylSeq; segTypeExt];
%     sylSeqTime = [sylSeqTime; timestamp];
%     sylTrans = [sylTrans; allSeq];
%     sylTransTime = [sylTransTime; allSeqTime'];
% end
% 
% %Step through the motif bands to calculate the recovery of sequencing
% for i = 1:size(bands,1)  
%     tempArray = sylTrans(bands(i,1):bands(i,2),1)';
%     
%     %Motif Fraction -- fraction of vocalizations that are in an identifyable motif
%     MF(i) = (length(strfind(tempArray, motif))*length(motif))/length(bands(i,1):bands(i,2));
%     
%     %Motif Continuity -- probability that once you start the sequence, you'll finish it
%     for j = 1:length(motif)
%         if j == 1
%             %This condition deals with variable repeats/stutters of the opening syllable
%             pntrs = strfind(tempArray, motif(1));
%             pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
%             MC(i, j) = length(pntrs_refined);
%         else
%         MC(i, j) = length(strfind(tempArray, motif(1:j)));
%         end
%     end
%     MC(i,:) = MC(i, :)./MC(i, 1);
% end
% 
% %It's some other day; keep array intact
% xs = [1, size(bands,1)];
% 
% %For summary stats
% mornBandI = bandTime < datenum([date mornEnd],'yymmddHH:MM:SS');
% nightBandI = bandTime > datenum([date nightStart],'yymmddHH:MM:SS');
% 
% h(5) = figure(5); clf
% %Plot Motif Fraction
% subplot(2,1,1); cla
% hold on
% for i = 1:size(xs,1)
%     plot(xs(i,1):xs(i,2), MF(xs(i,1):xs(i,2)), 'LineWidth', 3)
% end
% y = [0, 1];
% axis tight; ylim(y); 
% ylabel('Motif Fraction', 'FontSize', 10)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)
% 
% %Plot Motif Continuity
% subplot(2,1,2); cla
% hold on
% for i = 1:size(xs,1)
%     plot(xs(i,1):xs(i,2), MC(xs(i,1):xs(i,2), end), 'LineWidth', 3)
% end
% y = [0, 1];
% axis tight; ylim(y); 
% xlabel('Motif Renditions', 'FontSize', 10); ylabel('Motif Continuity', 'FontSize', 10)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)
% 
% set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7]);
% 
% %Calculate the interval mean and std of MF and MC for future stats
% meanMF = [mean(MF(mornBandI)), mean(MF(nightBandI))];
% stdMF = [std(MF(mornBandI)), std(MF(nightBandI))];
% 
% meanMC = [mean(MC(mornBandI,end)), mean(MC(nightBandI,end))];
% stdMC = [std(MC(mornBandI,end)), std(MC(nightBandI,end))];

%Save all data
save([outputDir, outputName, '_pharmaDataset.mat']);

%Save figures
%close all

display(['Completed ' annotName])

clear all
