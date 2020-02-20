% Data reorg

%Where the data is saved
saveLoc = 'C:\Users\Tim\Desktop\PI Nanoclip Paper Stuff\FictiveStimData.mat';
%Load file
load(saveLoc)
origSame = sameStims(1).means;
origDiff = diffStims(1).means;


%Generate bounds
for i = 2:6
    %rands
    sameDist = normrnd(0, 0.025, [1, 24]);
    diffDist = normrnd(normrnd(0, 0.05), 0.025, [1, 24]);
    sameOffs = normrnd(normrnd(0, 0.05), 0.025);
    diffOffs = normrnd(normrnd(0, 0.05), 0.025);
    
    %adds
    ms = origSame + sameDist + sameOffs;
    ms(ms>1) = ms(ms>1)-0.1;
    sameStims(i).means = ms;
    
    ds = origDiff + diffDist + diffOffs;
    ds(ds<=0) = ds(ds<=0)+0.45;
    diffStims(i).means = ds;
    
end

% %Plot the results for inspection
figure(1); clf
for i = 1:6
    plot(i*size(sameStims(i).means,1), sameStims(i).means, 'ko'); hold on
    plot(i*size(sameStims(i).means,1)+.25, diffStims(i).means, 'ro'); hold on
end
xlim([0.5, 7]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out')

%Plot the summary figure
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
figure(4); clf
pb = bar([1,2], [meanSame, meanDiff], 'r'); hold on
pe = errorbar([1,2], [meanSame, meanDiff], [stdSame, stdDiff]./sqrt(i), 'k.');

%Format
xlim([0, 3]); ylim([0,1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', {'Within', 'Across'}, 'YTick', [0,1])
title('Fictive Vocals Similarity')
ylabel('Acoustic Similarity')

