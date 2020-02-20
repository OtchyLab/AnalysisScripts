%Code written to perform the SAP-similarity analysis on the fictive
%singing recordings. Primary outputs are:
%   (1) some sort of low-D visualization of the data
%   (2) some defined scalar metric of the similarity/difference of vocalization
%
%This function will process an entire folder of *.mat files that each
%contain audio (and other) information about the vocalizations produced
%from a single set of stimulation parameters
%
%
%Written by TMO 01/30/20

%% Clear the workspace
clear

%Constants
fs = 44150;

%Source and destination data
% sLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\Fictive Syllable analysis'; %fictive data location
sLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\Naturally produced syllables'; %fictive data location

%Grab the file listing
fls = dir([sLoc, filesep, 'LR71*.mat']);
fnames = getFieldVectorCell(fls', 'name')';

%% Sequentially load and pre-process files
snips = [];
ID = [];
snipMS = 100; %in ms
snipSize = snipMS*fs/1000; %in samples
for i = 1:numel(fnames)
    %Assemble long name
    fn = [sLoc, filesep, fnames{i}];
    
    %Load from file
    S = load(fn);
    
    %Stretch params
    timeBins = size(S.data.template,2);
    rendLength = size(S.audio.aligned_audioCube,3);
    trimpnts = floor(S.data.templatesyllBreaks .* (fs/1000));
    maxfiles = min([numel(S.audio.raw), 30]);
    for j = 1:maxfiles
        
        %Create the linear warping path for the segment
        syllBreaks = getWarpedStarts(S.data.pFull{j},S.data.templatesyllBreaks);
        a = S.data.templatesyllBreaks'; path_t = [1; a(:); timeBins];
        b = syllBreaks; path_r = [1; b(:); rendLength];
        sp_path = [path_t,path_r];
        p = (1:timeBins)';
        q = interp1(path_t,path_r,1:timeBins,'linear')';
        
        %Do the warping for the audio timeseries
        try
            aligned_TS = alignRawTS(S.audio.raw{j}, rendLength, [p,round(q)]);
            
            %Trim to syllable boundaries
            snipCut = aligned_TS(trimpnts(1):trimpnts(2));
            
            %Linear warp to standardized length
            snipScaled = interp1(1:numel(snipCut), snipCut, linspace(1, numel(snipCut), snipSize));
            
            %Add it to the processing list
            snips = [snips; snipScaled];
            ID = [ID; i];
        catch
            display('Uh oh, skipped a file')
        end
    end
    
    %Delete the data loaded from file
    clear('S')
end

%% Run the pairwise similarity analysis
acc = [];
sim = [];
compID = [];
numSnips = size(snips,1);
for i = 1:numSnips
    i %to check in on where we be
    for j = 1:numSnips %cut comparisons in 1/2
        
        %Skip calculating similarity of identical snips
        if i ~= j
            [accScore, simScore] = SapSimilarity3(snips(i,:), snips(j,:), fs);
            if isnan(accScore)
                display('Uh... something is wrong...')
            end
            
            %Copy out the results into the carried array
            compID = [compID; ID(i), ID(j)];
            acc = [acc; accScore];
            sim = [sim; simScore];
            
        end
    end
end

%% Analyze the similarity results
sameM = []; sameS = [];
diffM = []; diffS = [];
stimTypes = sort(unique(compID(:,1)));

%Cycle trhough and process each syllable type
for i = 1:numel(stimTypes)
    %Setup the required masks
    sameStim = false(size(acc));
    diffStim = false(size(acc));
    
    curStim = stimTypes(i);
    
    %Same stim comparisons mask
    idx = find(compID(:,1)==curStim & compID(:,2)==curStim);
    sameStim(idx) = true;
    
    %Different stim comparisons mask
    idx = find(compID(:,1)==curStim & compID(:,2)~=curStim);
    diffStim(idx) = true;
    
    %Use the masks for descriptive stats
    sameM(i) = mean(acc(sameStim));
    sameS(i) = std(acc(sameStim));
    diffM(i) = mean(acc(diffStim));
    diffS(i) = std(acc(diffStim));

end

%% Save the output to commong file

%Where the data is going
% saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveStimData.mat'; %Fictive data location
saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveControlData.mat'; %Control data location
if exist(saveLoc)
    %Load file
    load(saveLoc)
    
    %Append new data to the end
    sameStims(end+1).means = sameM;
    %sameStims(end).stds = sameS;
    diffStims(end+1).means = diffM;
    %diffStims(end).stds = diffS;
    
else
    %Setup the variable to save
    sameStims.means = sameM;
    %sameStims.stds = sameS;
    diffStims.means = diffM;
    %diffStims.stds = diffS;
    
end

%Save data back to the carrying var on disk
save(saveLoc, 'sameStims', 'diffStims');


%% Load and plot the figures and stats for single dataset

%Load file
% saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveStimData.mat'; %Fictive data location
saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveControlData.mat'; %Control data location
load(saveLoc)

% %Plot the results for inspection
% figure(1); clf
% for i = 1:6
%     plot(i*size(sameStims(i).means,1), sameStims(i).means, 'ko'); hold on
%     plot(i*size(sameStims(i).means,1)+.25, diffStims(i).means, 'ro'); hold on
% end
% xlim([0.5, 7]); ylim([0,1])
% set(gca, 'Box', 'off', 'TickDir', 'out')

%Plot the paper figure
sameS = [];
diffS = [];

%Mean for each bird
for i = 1:6
    sameS = [sameS, mean(sameStims(i).means)];
    diffS = [diffS, mean(diffStims(i).means)];
end

%Mean across stimulation conditions
meanSame = mean(sameS); stdSame = std(sameS);
meanDiff = mean(diffS); stdDiff = std(diffS);

%Plot the figure
figure(2); clf
pb = bar([1,2], [meanSame, meanDiff], 'r'); hold on
pe = errorbar([1,2], [meanSame, meanDiff], [stdSame, stdDiff]./sqrt(i), 'k.');

%Format
xlim([0, 3]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', {'Within', 'Across'}, 'YTick', [0,1])
title('Fictive Vocals Similarity')
ylabel('Acoustic Similarity')

% Significance testing... paired t-test
[h1, p1] = kstest(sameS); % h1 == 1 (not normal)
[h2, p2] = kstest(diffS); % h2 == 1 (not normal)
[h3, p3] = ttest(sameS, diffS); % h = 1; p = 7.7869e-04 (normality assumption violated)
[p4, h4] = ranksum(sameS, diffS); % h = 1l p = 0.0043 (signif)

%% Load and plot the figures and stats for the fictive and control datasets

%Load file
saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveStimData.mat'; %Fictive data location
load(saveLoc)

%Collate the data for the fictive stim birds
sameStim = [];
diffStim = [];
for i = 1:6
    sameStim = [sameStim, mean(sameStims(i).means)]; %mean over syllables (within a bird)
    diffStim = [diffStim, mean(diffStims(i).means)]; %mean over syllables (within a bird)
end

%Mean across birds for each condition
meanSame = mean(sameStim); stdSame = std(sameStim);
meanDiff = mean(diffStim); stdDiff = std(diffStim);

%Load file
saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveControlData.mat'; %Control data location
load(saveLoc)

%Collate the data for the control birds
sameCntl = [];
diffCntl = [];
for i = 1:6
    sameCntl = [sameCntl, mean(sameStims(i).means)]; %mean over syllables (within a bird)
    diffCntl = [diffCntl, mean(diffStims(i).means)]; %mean over syllables (within a bird)
end

%Mean across birds for each condition
meanSameCntl = mean(sameCntl); stdSameCntl = std(sameCntl);
meanDiffCntl = mean(diffCntl); stdDiffCntl = std(diffCntl);

%Plot the figure
figure(3); clf
pb = bar(1:4, [meanSame, meanDiff, meanSameCntl, meanDiffCntl], 'b'); hold on
pe = errorbar(1:4, [meanSame, meanDiff, meanSameCntl, meanDiffCntl], [stdSame, stdDiff, stdSameCntl, stdDiffCntl]./sqrt(i), 'k.');

%Format
xlim([0, 5]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', {'Within', 'Across','Within', 'Across'}, 'YTick', [0,1])
title('Fictive Vocals Similarity')
xlabel('Fictive Vocals           Natural Vocals')
ylabel('Acoustic Similarity')

set(gcf, 'units', 'inches', 'position', [14.75, 8.25, 4.25, 4.5])

% Significance testing... paired t-test
[h1, p1] = kstest(sameStim); % h1 == 1 (not normal)
[h2, p2] = kstest(diffStim); % h2 == 1 (not normal)
[h3, p3] = kstest(sameCntl); % h3 == 1 (not normal)
[h4, p4] = kstest(diffCntl); % h4 == 1 (not normal)

% [h3, p3] = ttest(sameS, diffS); % h = 1; p = 7.7869e-04 (normality assumption violated)
[hFV1, pFV1] = ttest(sameStim, diffStim); % h = 1 p = 3.5e-5 (signif)
[hNV1, pNV1] = ttest(sameCntl, diffCntl); % h = 1 p = 0.0016 (signif)
[hCV1, pCV1] = ttest2(sameStim, sameCntl); % h = 1 p = 5.3e-5 (signif)
[hCV2, pCV2] = ttest2(diffStim, diffCntl); % h = 1 p = 8.4e-4 (signif)










