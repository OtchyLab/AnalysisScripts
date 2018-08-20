function processElectroFolder

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\');

%Get all data files that match the criterion
crit = '*electroDataset.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially process the datasets
mcSet25 = [];
for i = 1:length(filelist)
    dataName = filelist(i).name;
    
    %Fix the datasets
%       fixNumbering(dataName, dataFolder);
%       update4Counts(dataName, dataFolder);
%       includeMC(dataName, dataFolder);
%           includeMCnew(dataName, dataFolder)
          mcSet25 = getMCnew(dataName, dataFolder, mcSet25);
    
end
folder = 'C:\Users\Tim\Desktop\Nif Project Figures\Motif Continuity Analysis';
save([folder, filesep, 'Grn141 MCanal2.mat'], 'mcSet25')

display(['Finished with folder: ' dataFolder])


function mcSet = getMCnew(dataName, dataFolder, mcSet)
%Load data from file
load([dataFolder filesep dataName])
close all

mcSet = [mcSet; offset, finalMC];

function includeMCnew(dataName, dataFolder)
%Adds the motif continuity statistic to the pharamcologically lesioned datasets. Code is taken (and adapted) from the code
%for electrolytically lesioned data

%Load data from file
load([dataFolder filesep dataName])

%%%%%%%%%%%%%%%%%%%%%
% Recreate the transitions arrays without banding
%%%%%%%%%%%%%%%%%%%%%

%Reset Calculated Variables
bands2 = [];

%Calculate sequencing primitives for the entire annotation
[~, ~, sylTrans2, sylTransTime2] = simpleSeq(elements);

%Find the location and recording time of all motifs
motifLoc = strfind(sylTrans2(:,1)', motif);
motifTimes = sylTransTime2(motifLoc);

%For each time bin (morning, post-inj, and night), determine the indices of motifs that correspond to first/last 25 motifs
bands2(1,:) = [1, motifLoc(26)-1];                                                           %all syllables surrounding the 1st 25 motifs
bands2(3,:) = [motifLoc(end-26)+length(motif), length(sylTrans2)];             %all syllables surrounding the last 25 motifs 

if offset == 0
    %Find the index corresponding to the first syllable and motif post lesion
    lesionInx = find(sylTransTime2 > critTime,1, 'first');
    lesMotifInx = find((motifLoc-lesionInx) > 0, 1, 'first');
    
    bands2(2,:) = [lesionInx, motifLoc(lesMotifInx+26)-1];                          %all syllables surrounding the first 25 motifs post-lesion
end

%Step through the bands to calculate the recovery of sequencing
MC25 = [];
for i = 1:size(bands2,1)
    if bands2(i,1) ~= 0
        tempArray = sylTrans2(bands2(i,1):bands2(i,2),1)';
        
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
                MC25(i, j) = length(pntrs_refined);
            else
                MC25(i, j) = length(strfind(tempArray, motif(1:j)));
            end
        end
        MC25(i,:) = MC25(i, :)./MC25(i, 1);
    else
        MC25(i,:) = NaN(1,size(MC25,2));
    end
end

%Calculate the interval mean and std of MF and MC for future stats
finalMC = MC25(:,end)';

%%%%%%%%%%%%%%%%%%%%%
%Save data to file
%%%%%%%%%%%%%%%%%%%%%
%Save all data
save([dataFolder filesep dataName]);
% 
% %Save figures
% savefig(h, [outputDir, outputName, '.fig']);

display(['Completed ' annotName])

%Close any open windows
close all

function mcSet = getMC(dataName, dataFolder,mcSet)
%Load data from file
load([dataFolder filesep dataName])
close all

mcSet = [mcSet;offset, meanMC];

function includeMC(dataName, dataFolder2)
%Adds the motif continuity statistic to the pharamcologically lesioned datasets. Code is taken (and adapted) from the code
%for electrolytically lesioned data

%Load data from file
load([dataFolder2 filesep dataName])

% if exist('MC') ~=0 %check if present
%%%%%%%%%%%%%%%%%%%%%
%Set up sequencing and motif bands
%%%%%%%%%%%%%%%%%%%%%
% motif = [1,2,3];%Grn046
% motif = [2,3,4];%Grn121
motif = [2,3,4];%Grn141

date = char(sp{2});

mornEnd = '14:00:00';
nightStart = '19:00:00';

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
%Setup the windowing for analysis
% for j = 1:size(b,1)
%     %Calculate sequencing primitives
%     [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));
%     
% %     %Find the starts to all motifs
% %     Indx = strfind(allSeq(:,1)', motif);
% %     
% %     %Generate sequence start and ends
% %     try
% %         indWin = windowTS(Indx, winSize, 1, 'pad', 'boxcar');
% %     catch
% %         indWin = Indx;
% %     end
% %     motifBands = [indWin(:,1), indWin(:,end)+length(motif)-1];
% %     motifBands(isnan(motifBands(:,1)),1) = Indx(1);
% %     motifBands(isnan(motifBands(:,2)),2) = Indx(end)+length(motif)-1;
% %     motifBands(motifBands<1) = 1;
% %     motifBands(motifBands>size(allSeq,1)) = size(allSeq,1);
% %     bands = [bands;motifBands+size(sylTrans,1)];
% %     
% %     %Locate the center time for each band
% %     try
% %         t = allSeqTime(indWin(:,ceil(winSize/2)));
% %     catch
% %         t = allSeqTime(indWin(ceil(length(indWin)/2)));
% %     end
% %     bandTime = [bandTime; t'];
%     
%     %Push to permenant arrays
%     sylSeq = [sylSeq; segTypeExt];
%     sylSeqTime = [sylSeqTime; timestamp];
%     sylTrans = [sylTrans; allSeq];
%     sylTransTime = [sylTransTime; allSeqTime'];
% end

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

%Step through the motif bands to calculate the recovery of sequencing
%     for i = 1:size(bands,1)
%         tempArray = sylTrans(bands(i,1):bands(i,2),1)';
%
%         %Motif Fraction -- fraction of vocalizations that are in an identifyable motif
%         MF(i) = (length(strfind(tempArray, motif))*length(motif))/length(bands(i,1):bands(i,2));
%
%         %Motif Continuity -- probability that once you start the sequence, you'll finish it
%         for j = 1:length(motif)
%             if j == 1
%                 %This condition deals with variable repeats/stutters of the opening syllable
%                 pntrs = strfind(tempArray, motif(1));
%                 pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
%                 MC(i, j) = length(pntrs_refined);
%             else
%                 MC(i, j) = length(strfind(tempArray, motif(1:j)));
%             end
%         end
%         MC(i,:) = MC(i, :)./MC(i, 1);
%     end
%
%     %It's some other day; keep array intact
%     xs = [1, size(bands,1)];
%
%     %For summary stats
%     mornBandI = bandTime < datenum([date mornEnd],'yymmddHH:MM:SS');
%     nightBandI = bandTime > datenum([date nightStart],'yymmddHH:MM:SS');
%
%     h(5) = figure(5); clf
%     %Plot Motif Fraction
%     subplot(2,1,1); cla
%     hold on
%     for i = 1:size(xs,1)
%         plot(xs(i,1):xs(i,2), MF(xs(i,1):xs(i,2)), 'LineWidth', 3)
%     end
%     y = [0, 1];
%     axis tight; ylim(y);
%     ylabel('Motif Fraction', 'FontSize', 10)
%     set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)
%
%     %Plot Motif Continuity
%     subplot(2,1,2); cla
%     hold on
%     for i = 1:size(xs,1)
%         plot(xs(i,1):xs(i,2), MC(xs(i,1):xs(i,2), end), 'LineWidth', 3)
%     end
%     y = [0, 1];
%     axis tight; ylim(y);
%     xlabel('Motif Renditions', 'FontSize', 10); ylabel('Motif Continuity', 'FontSize', 10)
%     set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)
%
%     set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7]);

%Calculate the interval mean and std of MF and MC for future stats
%     meanMF = [mean(MF(mornBandI)), mean(MF(nightBandI))];
%     stdMF = [std(MF(mornBandI)), std(MF(nightBandI))];

meanMC = [mornMC(end), nightMC(end)];
%     stdMC = [std(MC(mornBandI,end)), std(MC(nightBandI,end))];

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

close all
%Do what needs to be updated
%%%%%%%%%%%%%%%\
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
xlim([0 400]); ylim([0,max([mornY, postY, nightY])*1.1])
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

%%%%%%%%%%%%%%%%%%%%%
%Set up sequencing and motif bands
%%%%%%%%%%%%%%%%%%%%%
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

%Update the critical timepoints for windowing (i.e., don't window across boundaries)
% b defines the start and end of chunks that won't be broken up/banded across
if offset == 0
    %Today is lesion day; mark the breakpoint
    X = find(time<critTime,1,'last');
    b = [1, X; X+1, numRends];
else
    %It's some other day; keep elements intact
    b = [1, numRends];
end

%Setup the windowing for analysis
for j = 1:size(b,1)
    %Calculate sequencing primitives
    [segTypeExt, timestamp, allSeq, allSeqTime] = simpleSeq(elements(b(j,1):b(j,2)));

    %Find the starts to all motifs
    Indx = strfind(allSeq(:,1)', motif);

    %Generate sequence start and ends
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

%Step through the motif bands to calculate the recovery of sequencing
for i = 1:size(bands,1)  
    tempArray = sylTrans(bands(i,1):bands(i,2),1)';
    
    %Motif Fraction -- fraction of vocalizations that are in an identifyable motif
    MF(i) = (length(strfind(tempArray, motif))*length(motif))/length(bands(i,1):bands(i,2));
    
    %Motif Continuity -- probability that once you start the sequence, you'll finish it
    for j = 1:length(motif)
        if j == 1
            %This condition deals with variable repeats/stutters of the opening syllable
            pntrs = strfind(tempArray, motif(1));
            pntrs_refined = pntrs(tempArray(pntrs+1) ~=  motif(1));
            MC(i, j) = length(pntrs_refined);
        else
        MC(i, j) = length(strfind(tempArray, motif(1:j)));
        end
    end
    MC(i,:) = MC(i, :)./MC(i, 1);
end

if offset == 0
    %Today is lesion day; mark the breakpoint
    xs = [1, find(bandTime < critTime,1,'last');...
           find(bandTime < critTime,1,'last')+1, size(bands,1)];
       
    %For summary stats
    mornBandI = (bandTime < critTime) & (bandTime < datenum([date mornEnd],'yymmddHH:MM:SS'));
    postBandI =  (bandTime > critTime) & (bandTime < (bandTime(find(bandTime>critTime,1,'first'))+(postDur/24))) & (bandTime < datenum([date nightStart],'yymmddHH:MM:SS'));
    nightBandI = (bandTime > critTime) & (bandTime > datenum([date nightStart],'yymmddHH:MM:SS'));
else
    %It's some other day; keep array intact
    xs = [1, size(bands,1)];
    
    %For summary stats
    mornBandI = bandTime < datenum([date mornEnd],'yymmddHH:MM:SS');
    postBandI = bandTime < -inf;
    nightBandI = bandTime > datenum([date nightStart],'yymmddHH:MM:SS');
end

h(4) = figure(4); clf
%Plot Motif Fraction
subplot(2,1,1); cla
hold on
for i = 1:size(xs,1)
    plot(xs(i,1):xs(i,2), MF(xs(i,1):xs(i,2)), 'LineWidth', 3)
end
y = [0, 1];
axis tight; ylim(y); 
ylabel('Motif Fraction', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

%Plot Motif Continuity
subplot(2,1,2); cla
hold on
for i = 1:size(xs,1)
    plot(xs(i,1):xs(i,2), MC(xs(i,1):xs(i,2), end), 'LineWidth', 3)
end
y = [0, 1];
axis tight; ylim(y); 
xlabel('Motif Renditions', 'FontSize', 10); ylabel('Motif Continuity', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches', 'Position', [0 0 6 7]);

%Calculate the interval mean and std of MF and MC for future stats
meanMF = [mean(MF(mornBandI)), mean(MF(postBandI)), mean(MF(nightBandI))];
stdMF = [std(MF(mornBandI)), std(MF(postBandI)), std(MF(nightBandI))];

meanMC = [mean(MC(mornBandI,end)), mean(MC(postBandI,end)), mean(MC(nightBandI,end))];
stdMC = [std(MC(mornBandI,end)), std(MC(postBandI,end)), std(MC(nightBandI,end))];

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