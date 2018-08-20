% function LongSongV2
%LongSongV2 Analyzes the behavior (i.e., song) for the ELECTROLYTICLY lesioned birds for the Nif project.
%It streamlines the analysis done in the original LongSong.m
% 3/29/2015

% Load data sources from file

%Nif Lesion Hits
% mother = 'V:\Grn046\';
% annotNames = {'Grn046_140515_annotation.mat';...
%                         'Grn046_140516_annotation.mat';...
%                         'Grn046_140517_annotation.mat';...
%                         'Grn046_140518_annotation.mat'};
% critTime = datenum(2014, 5, 15, 11, 32, 00); %5/15/2014 @ 11:32
% outputName = 'Grn046_test';
% motif = [1,2,3];
% cond = 'hit';

mother = 'V:\Grn121\';
annotNames = {'Grn121_140906_annotation.mat';...
                        'Grn121_140907_annotation.mat';...
                        'Grn121_140908_annotation.mat';...
                        'Grn121_140909_annotation.mat'};
critTime = datenum(2014, 9, 6, 14, 5, 00); %9/6/2014 2:05pm
outputName = 'Grn121_test';
motif = [2,3,4];
cond = 'hit';

% mother = 'V:\Grn141\';
% annotNames = {'Grn141_141210_annotation.mat';...
%                         'Grn141_141211_annotation.mat';...
%                         'Grn141_141212_annotation.mat';...
%                         'Grn141_141214_annotation.mat'};
% critTime = datenum(2014, 12, 10, 11, 55, 00); %12/10/14 11:55am
% outputName = 'Grn141_test';
% motif = [2,3,4];
% cond = 'hit';

% mother = 'V:\Sil450\';
% annotNames = {'Sil450_140906_annotation.mat';...
%                         'Grn121_140907_annotation.mat';...
%                         'Grn121_140908_annotation.mat';...
%                         'Grn121_140909_annotation.mat'};
% critTime = datenum(2014, 5, 15, 11, 32, 00); %5/15/2014 @ 11:32
% outputName = 'Sil450_test';
% motif = [1,2,3,4,5];

% 
% mother = 'V:\Grn010\';
% annotNames = {'Grn010_140324_annotationAuto.mat'};...
% critTime = datenum(2014, 3, 24, 23, 59, 59); %5/15/2014 @ 11:32
% outputName = 'Grn010_test';
% motif = [1,2,3,4,5,6];

% mother = 'V:\Grn089\';
% annotNames = {'Grn089_140629_annotation.mat'};...
% critTime = datenum(2014, 6, 29, 23, 59, 59); %5/15/2014 @ 11:32
% outputName = 'Grn089_test';
% motif = [3,4,2];

% mother = 'V:\Pur935\';
% annotNames = {'Grn121_140906_annotation.mat';...
%                         'Grn121_140907_annotation.mat';...
%                         'Grn121_140908_annotation.mat';...
%                         'Grn121_140909_annotation.mat'};
% critTime = datenum(2014, 6, 29, 23, 59, 59); %5/15/2014 @ 11:32
% outputName = 'Grn089_test';
% motif = [3,4,2];

% mother = 'V:\Grn011\';
% annotNames = {'Grn011_140313_annotation.mat';...
%                         'Grn011_140314_annotation.mat';...
%                         'Grn011_140315_annotation.mat';...
%                         'Grn011_140316_annotation.mat';...
%                         'Grn011_140317_annotation.mat'};
% critTime = datenum(2014, 3, 13, 10, 21, 00); %3/13/2014 @ 11:32
% outputName = 'Grn011_test';
% motif = [1,2,3,4];

%Define main vars
fs = 44150;
    %extracted
sylAbsStartTime = [];
sylStartTime = [];
sylEndTime = [];
sylType = [];
wavName = [];
daynights = [];

    %calculated
sylFeat = [];
sylTrans = [];
sylTransTime = [];
sylSeq = [];
sylSeqTime = [];
bands = [];
bandTime = [];
transMats = [];

%Process the annotations sequentially
numAnnots = length(annotNames);
for i = 1:numAnnots
    %Construct file folder name
    sp = regexp(annotNames(i), '_', 'split');
    date = char(sp{1}(2));
    folder = ['20' date(1:2) '-' date(3:4) '-' date(5:6)];
    fileSource = [mother folder filesep];
    
    %Get datenum for boundaries
    daynights(i,1) = datenum([date '00:00:00'],'yymmddHH:MM:SS');
    daynights(i,2) = datenum([date '23:59:59'],'yymmddHH:MM:SS');
    
    %Load the annotation from file
    elements = []; keys = [];
    load([mother annotNames{i}], 'elements', 'keys');
    numFiles = length(keys);
    
    %Cycle through renditions to retrieve relevant data
    for j = 1:numFiles
        sylAbsStartTime = [sylAbsStartTime; elements{j}.segAbsStartTimes'];
        sylStartTime = [sylStartTime; elements{j}.segFileStartTimes'];
        sylEndTime = [sylEndTime; elements{j}.segFileEndTimes'];
        sylType = [sylType; elements{j}.segType];
        
        %Capture the filename of all syllables
        numSyls = length(elements{j}.segType);
        ins = cell(1,numSyls);
        for k = 1:numSyls
            ins{k} = keys{j};
        end
        wavName = [wavName; ins'];

        %Unpack the acoustic features if desired
        bDoAcoustic = true;
        if (bDoAcoustic)
            %Calc the SAP features
            if strcmp(keys{j}(end), 'v')
                %It's a WAV
                features = koenigSpectral(audioread([fileSource keys{i}]), fs);
            else
                %It's a DAT
                [chans, ~] = getChannels([fileSource keys{i}]);
                features = koenigSpectral(chans(1,:), fs);
            end

            %Parse by syllables
            insE = zeros(1,numSyls);
            insP = zeros(1,numSyls);
            insAmp = zeros(1,numSyls);
            insAM = zeros(1,numSyls);
            insFM = zeros(1,numSyls);
            for k = 1:numSyls
                s = max([1, elements{j}.segFileStartTimes(k)*1000]);
                e = min([length(features.Entropy) elements{j}.segFileEndTimes(k)*1000]);
                snip = floor(s):ceil(e);
                insE(k) = mean(features.Entropy(snip));
                insP(k) = mean(features.Pitch(snip));
                insAmp(k) = mean(features.amplitude(snip));
                insAM(k) = mean(features.AM(snip));
                insFM(k) = mean(features.FM(snip));
            end
            sylFeat = [sylFeat; insE', insP', insAmp', insAM', insFM'];
        end
    end
    
    %Update the critical timepoints for windowing (i.e., don't window across boundaries)
    if (critTime > daynights(i,1)) && (critTime < daynights(i,2))
        %Crit point is this day; split the elements pre/post
        for j = 1:numFiles
            time(j) = getFileTime(keys{j});
        end
        X = find(time<critTime,1,'last');
        b = [1, X; X+1, numFiles];
    else
        %Crit point elsewhere; keep elements intact
        b = [1, numFiles];
    end
    

    %Setup the windowing for analysis
    for j = 1:size(b,1)
        %Calculate sequencing primitives
        [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));

        %Find the starts to all motifs
        Indx = strfind(allSeq(:,1)', motif);
        
        %Generate sequence start and ends
        winSize = 201;
        try
           indWin = windowTS(Indx, winSize, 1, 'pad', 'boxcar');
        catch
            indWin = Indx;
        end
        motifBands = [indWin(:,1), indWin(:,end)+length(motif)-1];
        motifBands(isnan(motifBands(:,1)),1) = Indx(1);
        motifBands(isnan(motifBands(:,2)),2) = Indx(end)+length(motif)-1;
        motifBands(motifBands<1) = 1;
        motifBands(motifBands>size(allSeq,1)) = size(allSeq,1);
        bands = [bands;motifBands+size(sylTrans,1)];

%         %Find the starts to all motifs
%         Indx = strfind(segTypeExt', motif);
% 
%         %Generate sequence start and ends
%         winSize = 51;
%         indWin = windowTS(Indx, winSize, 1, 'pad', 'boxcar');
%         motifBands = [indWin(:,1), indWin(:,end)+length(motif)-1];
%         motifBands(isnan(motifBands(:,1)),1) = Indx(1);
%         motifBands(isnan(motifBands(:,2)),2) = Indx(end)+length(motif)-1;
%         motifBands(motifBands<1) = 1;
%         motifBands(motifBands>length(segTypeExt)) = length(segTypeExt);
%         bands = [bands;motifBands+size(sylSeq,1)];

        %Locate the center time for each band
        try
            t = allSeqTime(indWin(:,ceil(winSize/2)));
        catch
            t = allSeqTime(indWin(ceil(length(indWin)/2)));
        end
        bandTime = [bandTime; t'];

        %Push to permenant arrays
        sylSeq = [sylSeq; segTypeExt];
        sylSeqTime = [sylSeqTime; timestamp];
        sylTrans = [sylTrans; allSeq];
        sylTransTime = [sylTransTime; allSeqTime'];
    
    end
    
    
end

%Calculate Syllable Durations
sylDur = (sylEndTime - sylStartTime) * 1000; %in ms

%Create syllable index
sylInx = (sylType>=1 & sylType<101) | sylType==103;

%%%%%%%%%%%%%%%%%%%%%
%Generate heatmaps for duration v entropy
%%%%%%%%%%%%%%%%%%%%%

if bDoAcoustic
    %Generate subsets
    subDur = sylDur(sylInx);
    subFeat = sylFeat(sylInx,1);
    subType = sylType(sylInx);
    subTime = sylAbsStartTime(sylInx);

    %Locate critcal breakpoints...
    pntr = find((critTime > daynights(:,1)) & (critTime < daynights(:,2))); %which day (of the series) does the lesion take place?
    critBreaks = [critTime;daynights(pntr:end,2)];

    %Use these breakpoints to locate the indices corresponding to those times.
   % critPntr = find(subTime<critTime,1,'last');
    for i = 1:length(critBreaks)
        critPntr(i) = find(subTime<critBreaks(i),1,'last');
    end
    
    %Generate the range of indices desires for each heatmap
    poolSize = 1000;
    snipInx = [];
    for i = 1:(length(critPntr)-1)
        snipInx = [snipInx; critPntr(i)-poolSize, critPntr(i)]; %Before the boundary
        snipInx = [snipInx; critPntr(i)+1, critPntr(i)+poolSize]; %After the boundary
    end
    snipInx = [snipInx; critPntr(end)-poolSize, critPntr(end)]; %Final boundary only has before data

    %Calculate and plot heatmaps for each of the subsets
    for i = 1:size(snipInx,1)
        subplot(size(snipInx,1),1,i)
        N(i,:,:) = ndhist(subFeat(snipInx(i,1):snipInx(i,2)), subDur(snipInx(i,1):snipInx(i,2)), 'edgesx', -4:0.1:0, 'edgesy', 10:5:400, 'prob', 'filter', 'axis', [-4 0 0 400]);
        if i == 1
            imagesc(squeeze(N(i,:,:))'); sc = caxis; colormap(jet); axis xy
        else
            imagesc(squeeze(N(i,:,:))', sc); colormap(jet); axis xy
        end
        set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [], 'LineWidth', 3, 'FontSize', 16)
     end
     set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [-4, -2, 0], 'XTick', [1:20:79], 'XTickLabels', [0, 100, 200, 300])
     set(gca, 'LineWidth', 3, 'FontSize', 16)
     xlabel('Syl Dur (ms)'); ylabel('log(Entropy)')
    
     subplot(size(snipInx,1),1,1)
     title([char(sp{1}(1)) ' D2 PDF'], 'FontSize', 16);
     
     set(gcf, 'Units', 'Inches');
     set(gcf, 'Position', [0 0 4 10])     

end

%Set the subset of syllable types to generate transition matrix
defaults = [100, 101, 103]; n = 1:10;
idx = ismember(n, sort(unique(sylSeq)));
sylTypes = sort([n(idx), defaults]);
numberSyllTypes=length(sylTypes);


%Step through the motif bands to calculate the recovery of sequencing
for i = 1:size(bands,1)
    %Generate transMatrix for each band
    transMatrix=zeros(numberSyllTypes,numberSyllTypes);   
    for j = bands(i,1):bands(i,2)
        from = find(sylTypes==sylTrans(j,1));
        to = find(sylTypes==sylTrans(j,2));
        transMatrix(from,to) = transMatrix(from,to)+1;
    end
    transMats(i, :, :) = transMatrix;

    %Sequence Linearity
    numUnits = length(unique(sylTrans(bands(i,1):bands(i,2), 1)));
    binMat(transMatrix>0) = 1;
    numTransTypes = sum(binMat);
    
    SL(i) = numUnits/numTransTypes;
    
    %Sequence Consistency
    [~, c] = ismember(motif,sylTypes);
    tab = [];
    for j = 1:(length(c)-1)
        tab(j) = transMatrix(c(j), c(j+1));
    end
    numMotifTrans = sum(tab);
    numTotalTrans = sum(transMatrix(:));
    
    SC(i) = numMotifTrans/numTotalTrans;
    
    %Motif Fraction
    MF(i) = (length(strfind(sylTrans(bands(i,1):bands(i,2),1)', motif))*length(motif))/length(bands(i,1):bands(i,2));
    
    %Motif Continuity
    for j = 1:length(motif)
        subCnt(i, j) = length(strfind(sylTrans(bands(i,1):bands(i,2),1)', motif(1:j)));
    end
    subCent(i,:) = subCnt(i, :)./subCnt(i, 1);

    %Motif 
    tempArray = sylTrans(bands(i,1):bands(i,2),1)';
    for j = 2:length(motif)
        indx = strfind(tempArray, motif(j));
        minusArray = tempArray(indx-1);
        seq(i,j-1) = length(find(minusArray == motif(j-1)))/length(indx);
    end
    
end

%Find index points for ploting boundaries on windowed data
critIdx = find(bandTime<critTime,1,'last');
for i = 1:size(daynights,1)
    nightIdx(i) = find(bandTime<daynights(i,2),1,'last');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h(1) = figure(1); clf
plotInd = [0, critIdx, nightIdx];

%Plot Sequence Consistency
subplot(3,1,1)
hold on
for i = 1:(length(plotInd)-1)
    plot([(plotInd(i)+1):(plotInd(i+1))], SC((plotInd(i)+1):(plotInd(i+1))), 'LineWidth', 3)
end
y = [0, 1];
line([critIdx, critIdx], y,'LineStyle', '--', 'Color', 'r')
axis tight; ylim([0,1]); 
ylabel('Seq Consistency', 'FontSize', 10)
% title(['Summary D1 Distance'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

%Plot Motif Fraction
subplot(3,1,2)
hold on
for i = 1:(length(plotInd)-1)
    plot([(plotInd(i)+1):(plotInd(i+1))], MF((plotInd(i)+1):(plotInd(i+1))), 'LineWidth', 3)
end
y = [0, 1];
line([critIdx, critIdx], y,'LineStyle', '--', 'Color', 'r')
axis tight; ylim([0,1]); 
ylabel('Motif Fraction', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

%Plot Motif Continuity
subplot(3,1,3)
hold on
for i = 1:(length(plotInd)-1)
    plot([(plotInd(i)+1):(plotInd(i+1))], subCent((plotInd(i)+1):(plotInd(i+1)), end), 'LineWidth', 3)
end
y = [0, 1];
line([critIdx, critIdx], y,'LineStyle', '--', 'Color', 'r')
axis tight; ylim([0,1]); 
ylabel('Motif Continuity', 'FontSize', 10)
xlabel('Rendition', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot heatmap of durations
%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter parameters
window = 7;
factor = 6;

%Other parameters
maxInterval = 400;
trialBinSize = 50;
intertapBinSize = 5;
upperScale = 7.5E-2;

TargetDur = sylDur;

intertapMatrix = zeros(length(TargetDur), maxInterval);
convMatrix = ones(trialBinSize, intertapBinSize);
for m = 1:length(TargetDur)
    if (TargetDur(m) > 0) && (TargetDur(m) < maxInterval)
        intertapMatrix(m,round(TargetDur(m))) = 1;
    end
end
matrix = conv2(intertapMatrix, convMatrix, 'same');

% Filter
hfilt = ones(window*factor, window)/(factor*window^2);
matrix = imfilter(matrix, hfilt, 'replicate');
    
% Normalize
for m = 1:size(matrix,1)
    matrix(m,:) = matrix(m,:)/sum(matrix(m,:));
end

lesionPnt = find(sylAbsStartTime<=critTime, 1, 'last');

h(2) = figure(2); clf
imagesc(matrix, [0 upperScale]);
hold on
line([0, 400], [lesionPnt, lesionPnt], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')

% Add day markers
for i = 1:(length(daynights(:,2))-1)
    index(i) = find(sylAbsStartTime >= daynights(i,2), 1, 'first');
    plot(-8, index(i), '>k', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'clip', 'off')
end

box off
set(gca, 'TickDir', 'out')
xlabel('Vocalization Duration (ms)', 'FontSize', 16)
ylabel('Vocalizations', 'FontSize', 16)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 6])
set(gca, 'LineWidth', 3, 'FontSize', 16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot duration distributions on the same axis
minval = min(sylDur);
maxval = max(sylDur);
bins = 100;
setSize = 500;

% Snips for recovery comparisons
snipIdx = [0, lesionPnt, index, length(sylDur)];
snipEpdfX = []; snipEpdfY = [];

%Pool all of the morning-of prelesion data together
snip = sylDur((snipIdx(1)+1):snipIdx(2));
[snipEpdfX(1,:), snipEpdfY(1,:)] = epdf(snip',bins,minval,maxval);

%Cycle through the other snips
for i = 2:(length(snipIdx)-1)
    %Select the subset of data for following the current boundary
    snip = sylDur((snipIdx(i)+1):min([snipIdx(i)+setSize, snipIdx(i+1)]));
    
    %Calculate epdf of the durations
    [snipEpdfX(end+1,:), snipEpdfY(end+1,:)] = epdf(snip',bins,minval,maxval);
    
    %Select the subset of data for preceding the next boundary
    snip = sylDur(max([snipIdx(i+1)-setSize, snipIdx(i)]):(snipIdx(i+1)));
    
    %Calculate epdf of the durations
    [snipEpdfX(end+1,:), snipEpdfY(end+1,:)] = epdf(snip',bins,minval,maxval);
    
end

h(3) = figure(3);
cla; hold on
for i = 1:size(snipEpdfX,1)
    patch([snipEpdfX(i,:), snipEpdfX(i,end), snipEpdfX(i,1)], [snipEpdfY(i,:) + (0.25*(i-1))+.25, (0.25*(i-1))+.25, (0.25*(i-1))+.25], 'r')
end
% axis tight
xlim([0 400]); ylim([.1, 2.55])
set(gca, 'TickDir', 'out', 'Box', 'off', 'YTick', '')
xlabel('Syl Dur (s)', 'FontSize', 16)
ylabel('Recovery -->', 'FontSize', 16)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 3 5])
set(gca, 'LineWidth', 2, 'FontSize', 16)
title('Recovery of Temporal Structure', 'FontSize', 12);


figure;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recovery of duration distribution over days
h(4) = figure(4);
cla; hold on
for i = 1:size(snipEpdfX,1)
    %Correlation with the pre-lesion distribution
    durCorr(i) = corr(snipEpdfY(1,:)', snipEpdfY(i,:)');
end
xs.all = 1:9;
xs.D = [1, 4, 6, 8]; % Day time
xs.N = [3, 5, 7, 9]; % Night time
xs.L = 2;
plot(xs.D, durCorr(xs.D), 'sb')
plot(xs.N, durCorr(xs.N), 'sk')
plot(xs.L, durCorr(xs.L), 'sr')
ylim([0.5, 1])
xlabel('Time', 'FontSize', 16)
ylabel('Correlation with Pre-Lesion', 'FontSize', 16)
set(gca, 'XTick', xs.all, 'XTickLabels', [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}])
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 4 3])
title('Recovery of Temporal Structure', 'FontSize', 12);

%%%%%%%%%%%%%%%%%%%%%%%
%Recovery of Sequencing over days
h(5) = figure(5); clf
setSize = 101;
xs.all = 1:9;
xs.D = [1, 4, 6, 8]; % Day time
xs.N = [3, 5, 7, 9]; % Night time
xs.L = 2;

%Sequence Continuity
subplot(3,1,1); cla
hold on
SCsnipsM = mean(SC((plotInd(1)+1):(plotInd(2))));
SCsnipsS = std(SC((plotInd(1)+1):(plotInd(2))));
for i = 2:(length(plotInd)-1)
    %Snip after current break
    snip = SC((plotInd(i)+1):min([(plotInd(i)+setSize), plotInd(i+1)]));
    SCsnipsM(end+1) = mean(snip);
    SCsnipsS(end+1) = std(snip);
    
    %Snip before next break
    snip = SC(max([(plotInd(i+1)-setSize), plotInd(i)]):plotInd(i+1));
    SCsnipsM(end+1) = mean(snip);
    SCsnipsS(end+1) = std(snip);
end
errorbar(xs.D, SCsnipsM(xs.D), SCsnipsS(xs.D), 'Color', 'b', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.N, SCsnipsM(xs.N), SCsnipsS(xs.N), 'Color', 'k', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.L, SCsnipsM(xs.L), SCsnipsS(xs.L), 'Color', 'r', 'Marker', 's', 'LineStyle', 'none')
%xlabel('Time', 'FontSize', 16)
ylabel('Seq Cons', 'FontSize', 16)
set(gca, 'XTick', xs.all , 'XTickLabels', [])
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 10)

%Motif Fraction
subplot(3,1,2); cla
hold on
MFsnipsM = mean(MF((plotInd(1)+1):(plotInd(2))));
MFsnipsS = std(MF((plotInd(1)+1):(plotInd(2))));
for i = 2:(length(plotInd)-1)
    %Snip after current break
    snip = MF((plotInd(i)+1):min([(plotInd(i)+setSize), plotInd(i+1)]));
    MFsnipsM(end+1) = mean(snip);
    MFsnipsS(end+1) = std(snip);
    
    %Snip before next break
    snip = MF(max([(plotInd(i+1)-setSize), plotInd(i)]):plotInd(i+1));
    MFsnipsM(end+1) = mean(snip);
    MFsnipsS(end+1) = std(snip);
end
errorbar(xs.D, MFsnipsM(xs.D), MFsnipsS(xs.D), 'Color', 'b', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.N, MFsnipsM(xs.N), MFsnipsS(xs.N), 'Color', 'k', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.L, MFsnipsM(xs.L), MFsnipsS(xs.L), 'Color', 'r', 'Marker', 's', 'LineStyle', 'none')
%xlabel('Time', 'FontSize', 16)
ylabel('Motif Frac', 'FontSize', 16)
set(gca, 'XTick', xs.all, 'XTickLabels', [])
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 10)

%Motif Continuity
subplot(3,1,3); cla
hold on
MCsnipsM = mean(subCent((plotInd(1)+1):(plotInd(2)), end));
MCsnipsS = std(subCent((plotInd(1)+1):(plotInd(2)), end));
for i = 2:(length(plotInd)-1)
    %Snip after current break
    snip = subCent((plotInd(i)+1):min([(plotInd(i)+setSize), plotInd(i+1)]), end);
    MCsnipsM(end+1) = mean(snip);
    MCsnipsS(end+1) = std(snip);
    
    %Snip before next break
    snip = subCent(max([(plotInd(i+1)-setSize), plotInd(i)]):plotInd(i+1), end);
    MCsnipsM(end+1) = mean(snip);
    MCsnipsS(end+1) = std(snip);
end
errorbar(xs.D, MCsnipsM(xs.D), MCsnipsS(xs.D), 'Color', 'b', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.N, MCsnipsM(xs.N), MCsnipsS(xs.N), 'Color', 'k', 'Marker', 's', 'LineStyle', 'none')
errorbar(xs.L, MCsnipsM(xs.L), MCsnipsS(xs.L), 'Color', 'r', 'Marker', 's', 'LineStyle', 'none')
xlabel('Time', 'FontSize', 16)
ylabel('Motif Comp', 'FontSize', 16)
set(gca, 'XTick', xs.all, 'XTickLabels', [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}])
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 5])

% %%%%%%%%%%%%%%%%%%%%%%
% %Save data to file
outputLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
outputFile = [char(sp{1}(1)) '_longitudinal'];
%Save all data
save([outputLocation, outputFile, '_allData V2.mat']);

%Save figures
savefig(h, [outputLocation, outputFile, '_figures V2.fig']);

%%%%%%%%%%%%%%%%%%%%%
% Save to output location
%%%%%%%%%%%%%%%%%%%%%

output = [];
output.name = char(sp{1}(1));
output.type = cond;
output.annots = annotNames;
output.xs = xs;
output.labels =  [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}];
output.durCorr = durCorr;
output.SCsnipsM = SCsnipsM;
output.SCsnipsS = SCsnipsS;
output.MFsnipsM = MFsnipsM;
output.MFsnipsS = MFsnipsS;
output.MCsnipsM = MCsnipsM;
output.MCsnipsS = MCsnipsS;

outPerm =  ['Longitudinal Summary V2.mat'];
m = exist([outputLocation, outPerm]);
if m ==2 
    %File already exists
    clear('longStats');
    load([outputLocation, outPerm], 'longStats')
    longStats(end+1) = output;
else
    %No file yet created
    longStats = output;
end

%Save the updated data to file
save([outputLocation, outPerm], 'longStats')

display('done')
