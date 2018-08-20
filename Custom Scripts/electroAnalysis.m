% function electroAnalysis
%4/8/2015
%This script is for analying 1D and 2D pdfs for song structure -- mainly for the ELECTROLYTIC
%Nif lesions and is generalized to collect summary info over many days.
%
%Update 05/15/15 to include higher-resolution heatmaps and all accompanying metrics
%Update 06/03/15 to include the new Motif Completion scripts that restricts analysis to 25 (or so) motifs.

clear all; close all

%%%%%%%%%%%%%%%%%%%%
%Set data sources and injection time
%%%%%%%%%%%%%%%%%%%%

annotName = 'Grn046_140522_annotation.mat';
critTime = datenum(2014, 5, 15, 11, 32, 00); %5/15/2014 @ 11:32
motif = [1,2,3];
cond = 'hit';

% annotName = 'Grn052_140617_annotation.mat';
% critTime = datenum(2014, 6, 14, 12, 56, 00); %6/14/14 @ 12:56pm
% motif = [2,3,4];
% cond = 'miss';

% annotName = 'Grn094_140619_annotation.mat';
% critTime = datenum(2014, 6, 19, 11, 44, 00); %6/19/14 @ 11:44am
% motif = [2,3,4,5];
% cond = 'miss';

% annotName = 'Grn121_140913_annotation.mat';
% critTime = datenum(2014, 9, 6, 14, 5, 00); %9/6/2014 2:05pm
% motif = [2,3,4];
% cond = 'hit';

% annotName = 'Grn141_141217_annotation.mat';
% critTime = datenum(2014, 12, 10, 11, 55, 00); %12/10/14 11:55am
% motif = [2,3,4];
% cond = 'hit';

% annotName = 'Grn165_150309_annotation.mat';
% critTime = datenum(2015, 3, 9, 12, 21, 00); %
% motif = [2,3,4];
% cond = 'miss';

% annotName = 'Grn186_150501_annotation.mat';
% critTime = datenum(2015, 4, 24, 11, 48, 00); %
% motif = [2,3,4,5];
% cond = 'hit';

%Generate all other file locations/names based on annotation name
mother = 'V:\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];
fileSource = [mother, char(sp{1}) filesep '20' char(sp{2}(1:2)) '-' char(sp{2}(3:4)) '-' char(sp{2}(5:6)) filesep];

%Load the annotation from file
load(annotLoc)

%Define output location
outputDir = ['C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\060315\' char(sp{1}) '\'];
outputName = [char(sp(1)) '_' char(sp(2))];

%Get datenum for boundaries
date = char(sp{2});

%%%%%%%%%%%%%%%%%%%%%
%Define main vars
%%%%%%%%%%%%%%%%%%%%%
fs = 44150;

    %extracted
sylAbsStartTime = [];
sylStartTime = [];
sylEndTime = [];
sylDur = [];
sylType = [];
wavName = [];

    %calculated
sylFeat = [];

%%%%%%%%%%%%%%%%%%%%%
%Strip Info from annotation
%%%%%%%%%%%%%%%%%%%%%
numRends = length(keys);
for i = 1:numRends
    %Unpack the time, sylType, and duration
    sylAbsStartTime = [sylAbsStartTime; elements{i}.segAbsStartTimes'];
    sylDur = [sylDur; [elements{i}.segFileEndTimes - elements{i}.segFileStartTimes]'];
    sylType = [sylType; elements{i}.segType];
    time(i) = getFileTime(keys{i});

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

delete('chans');
%%%%%%%%%%%%%%%%%%%%%
%Determine time relative to the lesion day
%%%%%%%%%%%%%%%%%%%%%
%Determine which day (relative to lesion) this is
sp = regexp(annotName, '_', 'split');
Day = sp{2};
dayDatenum = datenum(Day, 'yymmdd');
offset = dayDatenum - floor(critTime);

%%%%%%%%%%%%%%%%%%%%%
%Set up indices for epochs
%%%%%%%%%%%%%%%%%%%%%
%Create syllable indices
sylIdx = (sylType>=1 & sylType<101) | sylType==103;

mornEnd = '12:00:00';
nightStart = '19:00:00';

postDur = 1;
if offset == 0
    %This is lesion day
    mornIdx = (sylAbsStartTime < critTime) & (sylAbsStartTime < datenum([date mornEnd],'yymmddHH:MM:SS'));
    postIdx =  (sylAbsStartTime > critTime) & (sylAbsStartTime < (sylAbsStartTime(find(sylAbsStartTime>critTime,1,'first'))+(postDur/24))) & (sylAbsStartTime < datenum([date nightStart],'yymmddHH:MM:SS'));
    nightIdx = (sylAbsStartTime > critTime) & (sylAbsStartTime > datenum([date nightStart],'yymmddHH:MM:SS'));
else
    %Its some other day
    mornIdx = sylAbsStartTime < datenum([date mornEnd],'yymmddHH:MM:SS');
    postIdx = sylAbsStartTime < -inf;
    nightIdx = sylAbsStartTime > datenum([date nightStart],'yymmddHH:MM:SS');
end


%%%%%%%%%%%%%%%%%%%%%
%Calculate duration distributions
%%%%%%%%%%%%%%%%%%%%%
minval = 0;
maxval = 350;
binSize = 1.25; %bin size in ms
[mornX, mornY] = epdf_cbins(sylDur(sylIdx & mornIdx),binSize,minval,maxval);
[postX, postY] = epdf_cbins(sylDur(sylIdx & postIdx),binSize,minval,maxval);
[nightX, nightY] = epdf_cbins(sylDur(sylIdx & nightIdx),binSize,minval,maxval);
[allX, allY] = epdf_cbins(sylDur(sylIdx),binSize,minval,maxval);

%Plot the duration distributions on the same axis
h(1) = figure(1);
cla; hold on
plot(mornX,mornY, 'LineWidth', 2)
plot(postX,postY, 'LineWidth', 2)
plot(nightX,nightY, 'LineWidth', 2)

%Format the figure
axis tight
xlim([0 350]); ylim([0,max([mornY, postY, nightY])*1.1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 3, 'FontSize', 16)
xlabel('Syl Dur (ms)', 'FontSize', 16)
ylabel('P(t)', 'FontSize', 16)
set(gcf, 'Units', 'Inches', 'Position', [0 0 5 5])
title([sp{1} ' ' sp{2}], 'FontSize', 16);
legend({'morning'; 'post'; 'night'}); legend('boxoff')

%%%%%%%%%%%%%%%%%%%%%
%Generate heatmaps for duration v entropy
%%%%%%%%%%%%%%%%%%%%%
h(2) = figure(2); clf
if bDoAcoustic
    %Calculate and plot heatmaps for each of the subsets
    subplot(4,1,1); cla
    mornN = ndhist(sylFeat(sylIdx & mornIdx), sylDur(sylIdx & mornIdx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'axis', [-4 0 0 350]);
    imagesc(mornN'); sc = caxis; colormap(jet); axis xy; hold on
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [], 'LineWidth', 3, 'FontSize', 16)
    title([sp{1} ' ' sp{2} ' D2 PDF'], 'FontSize', 16);
    ylabel('Morning')
    
    if offset == 0
        subplot(4,1,2); cla
        postN = ndhist(sylFeat(sylIdx & postIdx), sylDur(sylIdx & postIdx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'axis', [-4 0 0 350]);
        imagesc(postN', sc); colormap(jet); axis xy; hold on
        set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [], 'LineWidth', 3, 'FontSize', 16)
        ylabel('Post')    
    end
    
    subplot(4,1,3); cla
    nightN = ndhist(sylFeat(sylIdx & nightIdx), sylDur(sylIdx & nightIdx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'axis', [-4 0 0 350]);
    imagesc(nightN', sc); colormap(jet); axis xy; hold on
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [], 'LineWidth', 3, 'FontSize', 16)
    ylabel('Night')
    
    subplot(4,1,4); cla
    allN = ndhist(sylFeat(sylIdx), sylDur(sylIdx), 'edgesx', -4:0.025:0, 'edgesy', 0:1.25:350, 'prob', 'axis', [-4 0 0 350]);
    imagesc(allN', sc); colormap(jet); axis xy; hold on
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [-4, -2, 0], 'XTick', [1:20:79], 'XTickLabels', [0, 100, 200, 300], 'LineWidth', 3, 'FontSize', 16)
    xlabel('Syl Dur (ms)'); ylabel('log(Entropy)')
    
     set(gcf, 'Units', 'Inches', 'Position', [0 0 4 10])     
end

%%%%%%%%%%%%%%%%%%%%%
%Generate heatmaps for duration over time
%%%%%%%%%%%%%%%%%%%%%
%The duration waterfall plot
h(3) = figure(3); clf
subplot(1,5,1:4); cla
[matrix, upperScale] = makeDurHeatmap(sylDur(sylIdx));
imagesc(matrix, [0 upperScale]);
y = ylim; colormap(jet);
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 1000:1000:size(matrix,1), 'LineWidth', 3, 'FontSize', 16)
xlabel('Syllable Durations (ms)', 'FontSize', 16)
title([sp{1} ' ' sp{2}], 'FontSize', 16);

%The colorbar showing snip locations
mornPatch = [find(mornIdx == true, 1, 'first'), find(mornIdx == true, 1, 'last')];
postPatch = [find(postIdx == true, 1, 'first'), find(postIdx == true, 1, 'last')];
nightPatch = [find(nightIdx == true, 1, 'first'), find(nightIdx == true, 1, 'last')];

subplot(1,5,5); cla
hold on
patch([0,0,1,1], [mornPatch fliplr(mornPatch)], [0, 0.45, 0.74], 'FaceAlpha', 1, 'EdgeColor', 'none')
if offset == 0
    patch([0,0,1,1], [postPatch fliplr(postPatch)], [0.85, 0.33, 0.1], 'FaceAlpha', 1, 'EdgeColor', 'none') 
end
patch([0,0,1,1], [nightPatch fliplr(nightPatch)], [0.93, 0.69, 0.13], 'FaceAlpha', 1, 'EdgeColor', 'none')
axis ij; axis tight
ylim(y);
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7]);

% %%%%%%%%%%%%%%%%%%%%%
% %Set up sequencing and motif bands
% %%%%%%%%%%%%%%%%%%%%%

%Reset Calculated Variables
bands = [];

%Motif # to caluclate over
numM = 25;

%Calculate sequencing primitives for the entire annotation
[~, ~, sylTrans, sylTransTime] = simpleSeq(elements);

%Find the location and recording time of all motifs
motifLoc = strfind(sylTrans(:,1)', motif);
motifTimes = sylTransTime(motifLoc);

%For each time bin (morning, post-inj, and night), determine the indices of motifs that correspond to first/last 25 motifs
bands(1,:) = [1, motifLoc(numM+1)-1];                                                           %all syllables surrounding the 1st numM motifs
bands(3,:) = [motifLoc(end-numM+1)+length(motif), length(sylTrans)];             %all syllables surrounding the last numM motifs 

if offset == 0
    %Find the index corresponding to the first syllable and motif post lesion
    lesionInx = find(sylTransTime > critTime,1, 'first');
    lesMotifInx = find((motifLoc-lesionInx) > 0, 1, 'first');
    
    bands(2,:) = [lesionInx, motifLoc(lesMotifInx+numM+1)-1];                          %all syllables surrounding the first numM motifs post-lesion
end

%Step through the bands to calculate the recovery of sequencing
MC = [];
for i = 1:size(bands,1)
    if bands(i,1) ~= 0
        tempArray = sylTrans(bands(i,1):bands(i,2),1)';
        
        %Motif Continuity -- probability that once you start the sequence, you'll finish it
        for j = 1:length(motif)
            if j == 1
                %This condition deals with variable repeats/stutters of the opening syllable
                pntrs = strfind(tempArray, motif(1));
                if pntrs(end) ~= length(tempArray)
                    pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
                else
                    pntrs_refined = pntrs(tempArray(pntrs(1:(end-1))+1) ~=  motif(1)); %deals with the unusual case in which the band ends with the opening syllable
                    pntrs_refined = [pntrs_refined, length(tempArray)];
                end
                MC(i, j) = length(pntrs_refined);
            else
                MC(i, j) = length(strfind(tempArray, motif(1:j)));
            end
        end
        MC(i,:) = MC(i, :)./MC(i, 1);
    else
        MC(i,:) = NaN(1,size(MC,2));
    end
end

%Calculate the interval mean and std of MF and MC for future stats
finalMC = MC(:,end)';


%%%%%%%%%%%%%%%%%%%%%
%Save data to file
%%%%%%%%%%%%%%%%%%%%%
%Save all data
save([outputDir, outputName, '_electroDataset.mat']);

%Save figures
savefig(h, [outputDir, outputName, '.fig']);

display(['Completed ' annotName])


%OLD CODE FROM PREVIOUS VERSIONS:
%
%This bit calculated the MC for large blocks of time; discarded 6/3/15
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
% %Update the critical timepoints for windowing (i.e., don't window across boundaries)
% % b defines the start and end of chunks that won't be broken up/banded across
% if offset == 0
%     %Today is lesion day; mark the breakpoint
%     X = find(time<critTime,1,'last');
%     b = [1, X; X+1, numRends];
% else
%     %It's some other day; keep elements intact
%     b = [1, numRends];
% end
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
% if offset == 0
%     %Today is lesion day; mark the breakpoint
%     xs = [1, find(bandTime < critTime,1,'last');...
%            find(bandTime < critTime,1,'last')+1, size(bands,1)];
%        
%     %For summary stats
%     mornBandI = (bandTime < critTime) & (bandTime < datenum([date mornEnd],'yymmddHH:MM:SS'));
%     postBandI =  (bandTime > critTime) & (bandTime < (bandTime(find(bandTime>critTime,1,'first'))+(postDur/24))) & (bandTime < datenum([date nightStart],'yymmddHH:MM:SS'));
%     nightBandI = (bandTime > critTime) & (bandTime > datenum([date nightStart],'yymmddHH:MM:SS'));
% else
%     %It's some other day; keep array intact
%     xs = [1, size(bands,1)];
%     
%     %For summary stats
%     mornBandI = bandTime < datenum([date mornEnd],'yymmddHH:MM:SS');
%     postBandI = bandTime < -inf;
%     nightBandI = bandTime > datenum([date nightStart],'yymmddHH:MM:SS');
% end
% 
% h(4) = figure(4); clf
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
% meanMF = [mean(MF(mornBandI)), mean(MF(postBandI)), mean(MF(nightBandI))];
% stdMF = [std(MF(mornBandI)), std(MF(postBandI)), std(MF(nightBandI))];
% 
% meanMC = [mean(MC(mornBandI,end)), mean(MC(postBandI,end)), mean(MC(nightBandI,end))];
% stdMC = [std(MC(mornBandI,end)), std(MC(postBandI,end)), std(MC(nightBandI,end))];

