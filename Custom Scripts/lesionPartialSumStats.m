function [sumFrac, corrReduct, powReduct] = lesionPartialSumStats(partialSumStats)

%Takes the output of LongSong.m
summary = [];

%Extract from the carried structure
name = getStructField(partialSumStats,'name');
date = getStructField(partialSumStats,'date');
chan = getStructField(partialSumStats,'chan');
frac = getStructField(partialSumStats,'lesionFrac');
powMeans = getStructField(partialSumStats,'powMeans');
recoveryMeans = getStructField(partialSumStats,'corrRecovery');

%Gen the list of all birds
birds = unique(name);

%Cycle through each bird and extract "best" example and add summary data to the running list
for i = 1:size(birds,2);
    %Gen mask for individual bird
    birdMask = strcmp(birds(i), name);
    
    subPow = powMeans(birdMask,:);
    subDate = date(birdMask);
    subRec = recoveryMeans(birdMask,:);
    subChan = chan(birdMask);
    subFrac = frac(birdMask);
    
    %Locate the best/most dramatic correlation drop
    normRec = [];
    for j = 1:size(subRec,1)
        normRec(j,:) = subRec(j,:)./subRec(j,1);
    end
    [minCor, idxBest] = min(normRec(:,2));
    summary(end+1).name = birds(i);
    summary(end).date = subDate(idxBest);
    summary(end).chan = subChan(idxBest);
    summary(end).pow = subPow(idxBest,:);
    summary(end).rec = subRec(idxBest,:);
    summary(end).frac = subFrac(idxBest);
end

%Plot summary figure
sumName = getStructField(summary,'name');
sumPow = getStructField(summary,'pow');
sumRec = getStructField(summary,'rec');
sumFrac = getStructField(summary,'frac');

%Normalize recovery by pre-lesion
for i = 1:size(sumRec,1)
    normSumPow(i,:) = sumPow(i,:)./sumPow(i,1);
%     normSumRec(i,:) = sumRec(i,:)./sumRec(i,1);
    normSumRec(i,:) = sumRec(i,:);
end
corrReduct = normSumRec(:,2);
powReduct = normSumPow(:,2);

%Best fit lines
[corrFit, corrStats] = polyfit(sumFrac,corrReduct,1);
[powFit, powStats] = polyfit(sumFrac,powReduct,1);
xss = [1,100];

%Plot it
figure(1);clf
subplot(2,2,1); cla
for i = 1:size(sumRec,1)
    plot([1:3],sumRec(i,:),'--o', 'DisplayName', birds{i}); hold on
end
xlim([0.5, 3.5]); ylim([0, 1])
ylabel('Mean Correlation to Pre')
set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', 1:3, 'XTickLabel', [{'Pre'}, {'Post'}, {'Night'}], 'YTick', 0:.5:1)

subplot(2,2,2); cla
scatter(sumFrac,corrReduct','o'); hold on
plot(xss, polyval(corrFit, xss), '--k')
xlim([0, 100]); ylim([0, 1])
axis square
xlabel('Lesion Size'); ylabel('Correlation')
set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', 0:50:100, 'YTick', 0:.5:1)

subplot(2,2,3); cla
for i = 1:size(sumPow,1)
    plot([1:3],sumPow(i,:),'--o', 'DisplayName', birds{i}); hold on
end
xlim([0.5, 3.5]); ylim([0, 0.25])
ylabel('HVC Activity (Norm to Pre)')
set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', 1:3, 'XTickLabel', [{'Pre'}, {'Post'}, {'Night'}], 'YTick', 0:.5:1)

subplot(2,2,4); cla
scatter(sumFrac,powReduct','o'); hold on
plot(xss, polyval(powFit, xss), '--k')
xlim([0, 100]); ylim([0, 1])
axis square
xlabel('Lesion Size'); ylabel('Reduction in HVC Activity (%)')
set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', 0:50:100, 'YTick', 0:.5:1)


%Significance of trend testing
[r.corr, p.corr] = corr(sumFrac, corrReduct); % r = -0.9114; p = 9.4292e-05
[r.pow, p.pow] = corr(sumFrac, powReduct)  % r = -0.8719; p = 4.6765e-04





