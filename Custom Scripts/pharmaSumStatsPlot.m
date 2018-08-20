function pharmaSumStatsPlot(pharmaStats)

summary = [];

%Grab info
name = getStructField(pharmaStats,'name');
date = getStructField(pharmaStats,'date');
offset = getStructField(pharmaStats,'offset');
d2emd = getStructField(pharmaStats,'D2emd');

birds = unique(name);
conditions = pharmaStats(1).conditions;

%Assemble figures for individuall birds
for i = 1:length(birds)
    birdMask = strcmp(birds(i), name);
    
    %X-vals for plotting
    notMT = [];
    for j = 1:length(offset{birdMask})
        notMT(end+1) = ~isempty(offset{birdMask}{j});
    end
    notMT = logical(notMT);
    notMT = notMT([1:3, 7:11]); % Strip down for use in subSets
    
    xs = [1:3, 7:11]-6;
    xs = xs(notMT);
    x_lim = [-6, 6];
    
    %Bar plot all of the similarity/distance measures on a single figure
    h(i) = figure(i); clf

    %Bar plot the 2D emd
    subSet = d2emd(birdMask,notMT);
    hold on
    plot(xs, subSet, 'Color', 'k', 'LineStyle', '--', 'Marker', 'o');
         %format
    xlim(x_lim); %axis square
    ylabel('Distance', 'FontSize', 8)
    title([char(birds(i)) ' D2 EMD'], 'FontSize', 8);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [-6:2:6], 'LineWidth', 2, 'FontSize', 10)
     
end

%Plot mean features son a single figure
h(end+1) = figure(20); clf
x_lims = size(subSet,2);
x_lim = [0,5];

%Bar plot the 2D emd
subSet = real(subIt(d2emd, name)); subSet(subSet==0) = NaN;
x_lims = size(subSet,2);
% subplot(2,7,10)
hold on
errorbar(1:x_lims, nanmean(subSet,1), nanstd(subSet,1), 'Color', 'k', 'LineStyle', '--', 'Marker', 'o');
for i = 1:size(subSet,1)
    plot(1:x_lims, subSet(i,:), 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);
end
plot([0,0], ylim,':r')
     %format
xlim(x_lim); %axis square
ylabel('Distance', 'FontSize', 8)
title([' D2 EMD'], 'FontSize', 8);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [ {'L-5d'}, {'1st'}, {'2nd'}, {'L+5d'}], 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches', 'Position', [0 0 6 3])


%Grab subsets for this treatment
subdist = subIt(d2emd, name); subdist(subdist==0) = NaN; %this is actually the EMD!!!!

%Normalize to baseline variation/difference
% for i = 1:size(subKL,1)
%     subKL(i,:) = subKL(i,:)./subKL(i,1);
%     subdist(i,:) = subdist(i,:)./subdist(i,1);
% end
  
%Plot polar plot of Summary differences in 2 metrics
h(end+1) = figure(102); clf
d(1:4) = plot(ones(size(subdist,1),1)*[1:4], subdist, '.k'); hold on
eb(2) = errorbar(1:4, mean(subdist,1), std(subdist,1)./sqrt(size(subdist,1)),'sk');
xlim([0.5,4.5]); ylim([0,40]); %axis square
ylabel('EMD', 'FontSize', 12)
title(['Earth Mover Distance'], 'FontSize', 12);
set(d, 'Marker', '.', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'LineStyle', 'none');
set(eb(2), 'Marker', 's', 'Color', 'k', 'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'L-5d'}, {'1st'}, {'2nd'}, {'L+5d'}], 'YTick', [0:5], 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches', 'Position', [0 0 6 3])

%Hypothesis Testing (1 day song diff from 5 day song)
[testCntrl.d1,pCntrl.d1,ciCntrl.d1,~] = ttest2([subdist(:,1)],[subdist(:,2)]);
[testCntrl.d2,pCntrl.d2,ciCntrl.d2,~] = ttest2([subdist(:,1)],[subdist(:,3)]);
[testCntrl.d5,pCntrl.d5,ciCntrl.d5,~] = ttest2([subdist(:,1)],[subdist(:,4)]);

%Save figure and data
saveFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Pharma\';
saveName = 'Final Pharma Stats 05152015';

%Save figures
savefig(h, [saveFolder, saveName '.fig']);
close all; clear('h')

%Save summary data
save([saveFolder, saveName '.mat']);

function subSet = subIt(feat, name)
%Do nothing; simply pass through
%subSet = feat;

%Parse data for 1st, 2nd, last post.
subSet = [];
for i = 1:size(feat,1)
    if strcmp(name(i),'Grn010')
        subSet(i,:) = feat(i, [1, 7, 10, 11]);
    elseif strcmp(name(i),'Grn011')
        subSet(i,:) = feat(i, [1, 8, 9, 11]);
    elseif strcmp(name(i),'Grn089')
        subSet(i,:) = feat(i, [1, 8, 9, 11]);
    elseif strcmp(name(i),'Grn091')
        subSet(i,:) = feat(i, [1, 8, 9, 11]);
    elseif strcmp(name(i),'Pur935')
        subSet(i,:) = feat(i, [1, 8, 9, 11]);
%     elseif strcmp(name(i),'Pur949') %Remove this line to exclude the last bird from summary stats.
%         subSet(i,:) = feat(i, [7, 8, 9]-4);
    end
    
end


