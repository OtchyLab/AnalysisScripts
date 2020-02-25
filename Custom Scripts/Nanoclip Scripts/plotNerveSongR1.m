%Script plots the acoustic similarity data charting nerve function recovery
%for the nanoclip paper. Also computes the required significance tests as
%reported in the paper. These figures/this analysis is in response to
%reviewer requests from NatComms.
%
%The script reads in the rawdata from a previously
%constructed .mat file containing all of the data from each condition.
%
% This script created the image submitted in the revised ms.
%
% Written by TMO 02-11-20

%Clear it all
clear all

%Load from file
load('/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/nerveBirdsData.mat') % <-- this is the newest data
if ~exist('summaryData', 'var')
    disp('Uh-oh... missing data?')
    return
end

%% Plot the individual bird traces
figure(100); clf
set(gcf,'Units', 'Inches', 'Position', [6, 5.25, 7.75, 6])

%Plot constants
col = {'b'; 'k'; 'r'; 'g'}; %{Nanoclip, Sham, Crush, Intact}
days = [-3:-1,1:8];

%Plot clip birds
block = 1:3;
b = [];all_b = [];
for i = block
    plot(days, summaryData(i).MeanSpect, ':b', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b];
mS(1) = mean(mean(b(1:3,:), 2));
means(1,:) = mean(b, 2)'./mS(1);
stds(1,:) = std(b, 1, 2)'./mS(1);

%Plot isolation birds
block = 4:6;
b = [];
for i = block
    plot(days, summaryData(i).MeanSpect, ':.k', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b];
mS(2) = mean(mean(b(1:3,:), 2));
means(2,:) = mean(b, 2)'./mS(2);
stds(2,:) = std(b, 1, 2)'./mS(2);

%Plot crush birds
block = 7:9;
b = [];
for i = block
     plot(days, summaryData(i).MeanSpect, ':.r', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b];
mS(3) = mean(mean(b(1:3,:), 2));
means(3,:) = mean(b, 2)'./mS(3);
stds(3,:) = std(b, 1, 2)'./mS(3);

%Plot intact birds
block = 10:12;
b = [];
for i = block
     plot(days, summaryData(i).MeanSpect, ':.g', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b];
mS(4) = mean(mean(b(1:3,:), 2));
means(4,:) = mean(b, 2)'./mS(4);
stds(4,:) = std(b, 1, 2)'./mS(4);

%Format figure
xlim([-3.5, 8.5]); ylim([0.5, 1.05])
xlabel('Time (days)'); ylabel('Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0.6:0.1:1)

%% Plot the timeseries data traces
figure(101); clf
set(gcf, 'Units', 'Inches', 'Position', [6, 6.25, 6.25, 5])   
pre = 1:3; post = 4:11;
for i=1:4
    %Plot pre
    jit = normrnd(0, 0.05, [1, numel(pre)]);
%     shadedErrorBar(days(pre), means(i, pre), stds(i, pre), col{i}, 1); hold on
    errorbar(days(pre)+jit, means(i, pre), stds(i, pre), col{i}); hold on
    
    %Plot post
    jit = normrnd(0, 0.05, [1, numel(post)]);
%     shadedErrorBar(days(post), means(i, post), stds(i, post), col{i}, 1);
    errorbar(days(post)+jit, means(i, post), stds(i, post), col{i});
end

%Format figure
xlim([-3.5, 8.5]); ylim([0.75, 1.05])
xlabel('Time (days)'); ylabel('Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -2:2:8, 'YTick', 0.75:0.25:1)

%% Plot the similarity differences
figure(102); clf
set(gcf, 'Units', 'Inches', 'Position', [6, 6.25, 6.25, 5])
base = 1:3; d1 = 4; d7 = 11;

%Plot each experimental group
deltaM = []; deltaS = [];
for i = 1:4 %3 groups
    ds = [];
    for j = 1:3 %
        pntr = 3*(i-1)+j;
        curBase = mean(all_b(1:3, pntr));
        ds = [ds; curBase-all_b(d1, pntr), curBase-all_b(d7, pntr)];
    end
    deltaM = [deltaM, ds(:,1), ds(:,2)];
%     deltaS = [deltaS, std(ds(:,1)), std(ds(:,2))];
end

bar(1:size(deltaM,2), mean(deltaM,1)); hold on
errorbar(1:size(deltaM,2), mean(deltaM,1), std(deltaM,1), '.k');

%Format figure
xlim([0, 9]); ylim([-0.025, 0.25])
xlabel('Groups'); ylabel('\Delta Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:8, 'YTick', [0, 0.25])


%% Compute significance tests for the time series

%Group Indices
ind = [1:3; 4:6; 7:9; 10:12];
tp = 1:numel(days);
nidx = [1,1,1,2,2,2,3,3,3,4,4,4];

%Normalize each bird data
rawdata = cell2mat(getFieldVectorCell(summaryData, 'MeanSpect'))';
for i = 1:numel(nidx)
    %normalize by group means
    rawdata(i,:) = rawdata(i,:)./mS(nidx(i));
end

%Check that all group are from normal distributions (for t-tests)
ksStats = [];
for i = tp
    for j = 1:4
        ksStats(j,i) = kstest(rawdata(ind(j,:),tp));
    end
end
allNormal = all(ksStats(:)); % allNormal = true

%Perform paired t-tests of the mean-pre to every subsequent day of that group
within_hSignif = []; within_pSignif = [];
for i = tp
    for j = 1:4
         base = mean(rawdata(ind(j,:),1:3),2);
        [within_hSignif(j,i), within_pSignif(j,i)] = ttest(base, rawdata(ind(j,:),i));
    end
end

%Perform unpaired t-tests of sham group to exp groups on every subsequent day of that group
across_hSignif = []; across_pSignif = [];
for i = tp
    s = 1;
    for j = [1,2,3]
        base = rawdata(ind(4,:),i);
        [across_hSignif(s,i), across_pSignif(s,i)] = ttest2(base, rawdata(ind(j,:),i));
        s = s+1;
    end
end
a = 1;
%Hardcoded summary:
%
% allNormal = true
% 
% within_hSignif =      0     0     0     0     0     0     0     0     1     0     0
%                       0     0     0     0     0     0     0     0     0     0     0
%                       0     0     0     1     1     1     1     1     1     1     1
%                       0     0     0     0     0     0     0     0     0     0     0
% 
% within_pSignif =      0.1356    0.1580    0.0580    0.0730    0.0668    0.0652    0.1512    0.3613   0.0390    0.0752    0.4175
%                       0.1465    0.5990    0.1004    0.5332    0.7207    0.4620    0.4719    0.4695   0.1939    0.4983    0.2786
%                       0.9081    0.5092    0.4649    0.0204    0.0291    0.0262    0.0223    0.0384   0.0404    0.0400    0.0166
%                       0.0588    0.2632    0.0725    0.7001    0.0805    0.0545    0.1096    0.4409   0.2609    0.1826    0.7315
%     
% across_hSignif =      0     0     0     0     0     0     0     0     0     0     0
%                       0     0     0     0     0     0     0     0     0     0     0
%                       0     0     0     1     1     1     1     1     1     1     1
%     
% 
% across_pSignif =      0.7137    0.9236    0.7335    0.1784    0.2008    0.2358    0.4766    0.6993    0.2709    0.2615    0.5126
%                       0.9367    0.7530    0.6176    0.4863    0.5617    0.5352    0.4688    0.7060    0.5979    0.9084    0.5167
%                       0.7349    0.9684    0.6038    0.0019    0.0019    0.0013    0.0006    0.0069    0.0073    0.0049    0.0088

%% Compute significance tests for the summary data

%Test for normality
ksTest2 = [];
for i = 1:size(deltaM,2)
    ksStats2(i) = kstest(deltaM(:,i));
end
allNormal2 = all(ksStats2(:)); % allNormal = false

%Test for difference from mean == 0
hmean0 = []; pmean0 = [];
for i = 1:size(deltaM,2)
    [hmean0(i), pmean0(i)] = ttest(deltaM(:,i));
end
a = 1;

%Hardcoded summary:
%
% allNormal2 = true
% 
% hmean0 =      0     0     0     0     1     1     0     0
%
% pmean0 =      0.0730    0.4175    0.5332    0.2786    0.0204    0.0166    0.7002    0.7315
