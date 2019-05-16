%Script plots the acoustic similarity data charting nerve function recovery
%for the nanoclip paper. Also computes the required significance tests as
%reported in the paper. The script reads in the rawdata from a previously
%constructed .mat file containing all of the data from each condition.
%
% Written by TMO 05-14-19

%Clear it all
clear all

%Load from file
load('/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/nerveBirdsData.mat')
if ~exist('summaryData', 'var')
    disp('Uh-oh... missing data?')
    return
end

%% Plot the individual bird traces
figure(100); clf
set(gcf,'Units', 'Inches', 'Position', [6, 5.25, 7.75, 6])

%Plot constants
col = {'b'; 'k'; 'r'}; %{Nanoclip, Sham, Crush}
days = [-3:-1,1:8];

%Plot clip birds
block = 1:3;
b = [];all_b = [];
for i = block
    plot(days, summaryData(i).MeanSpect, ':b', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b'];
mS(1) = mean(mean(b(1:3,:), 2));
means(1,:) = mean(b, 2)'./mS(1);
stds(1,:) = std(b, 1, 2)'./mS(1);
% shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'b', 1); hold on

%Plot isolation birds
block = 4:6;
b = [];
for i = block
    plot(days, summaryData(i).MeanSpect, ':.k', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b'];
mS(2) = mean(mean(b(1:3,:), 2));
means(2,:) = mean(b, 2)'./mS(2);
stds(2,:) = std(b, 1, 2)'./mS(2);
% shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'k', 1)

%Plot crush birds
block = 7:9;
b = [];
for i = block
     plot(days, summaryData(i).MeanSpect, ':.r', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
all_b = [all_b, b'];
mS(3) = mean(mean(b(1:3,:), 2));
means(3,:) = mean(b, 2)'./mS(3);
stds(3,:) = std(b, 1, 2)'./mS(3);
% shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'r', 1)

%Format figure
xlim([-3.5, 8.5]); ylim([0.5, 1.05])
xlabel('Time (days)'); ylabel('Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0.6:0.1:1)



%% Plot the summary data traces
figure(101); clf
set(gcf,'Units', 'Inches', 'Position', [6, 5.25, 6.25, 6])    
pre = 1:3; post = 4:11;
for i=1:size(all_b,1)
    %Plot pre
    shadedErrorBar(days(pre), means(i, pre), stds(i, pre), col{i}, 1); hold on
    
    %Plot post
    shadedErrorBar(days(post), means(i, post), stds(i, post), col{i}, 1);
end

%Format figure
xlim([-3.5, 8.5]); ylim([0.7, 1.05])
xlabel('Time (days)'); ylabel('Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -2:2:8, 'YTick', 0.75:0.25:1)



%% Compute significance tests

%Group Indices
ind = [1:3; 4:6; 7:9];
tp = 1:numel(days);
nidx = [1,1,1,2,2,2,3,3,3];

%Normalize each bird data
rawdata = cell2mat(getFieldVectorCell(summaryData, 'MeanSpect'))';
for i = 1:numel(nidx)
    %normalize by group means
    rawdata(i,:) = rawdata(i,:)./mS(nidx(i));
end

%Check that all group are from normal distributions (for t-tests)
ksStats = [];
for i = tp
    for j = 1:3
        ksStats(j,i) = kstest(rawdata(ind(j,:),tp));
    end
end
allNormal = all(ksStats(:)); % allNormal = true

%Perform paired t-tests of the mean-pre to every subsequent day of that group
within_hSignif = []; within_pSignif = [];
for i = tp
    for j = 1:3
         base = mean(rawdata(ind(j,:),1:3),2);
        [within_hSignif(j,i), within_pSignif(j,i)] = ttest(base, rawdata(ind(j,:),i));
    end
end

%Perform unpaired t-tests of sham group to exp groups on every subsequent day of that group
across_hSignif = []; across_pSignif = [];
for i = tp
    s = 1;
    for j = [1,3]
        base = rawdata(ind(2,:),i);
        [across_hSignif(s,i), across_pSignif(s,i)] = ttest2(base, rawdata(ind(j,:),i));
        s = s+1;
    end
end

%Hardcoded summary:
%
% allNormal = true
% 
% within_hSignif =      0     0     0     0     0     0     0     0     1     0     0
%                       0     0     0     0     0     0     0     0     0     0     0
%                       0     0     0     1     1     1     1     1     1     1     1
% 
% within_pSignif =      0.1356    0.1580    0.0580    0.0730    0.0668    0.0652    0.1512    0.3613   0.0390    0.0752    0.4175
%                       0.1465    0.5990    0.1004    0.5332    0.7207    0.4620    0.4719    0.4695   0.1939    0.4983    0.2786
%                       0.9081    0.5092    0.4649    0.0204    0.0291    0.0262    0.0223    0.0384   0.0404    0.0400    0.0166
%     
% across_hSignif =      0     0     0     0     0     0     0     0     0     0     0
%                       0     0     0     1     1     1     1     1     1     1     1
%     
% 
% across_pSignif =      0.7556    0.5975    0.8655    0.3393    0.2277    0.4146    0.3345    0.8786   0.3160    0.3911    0.9726
%                       0.7716    0.7393    0.8977    0.0026    0.0014    0.0014    0.0003    0.0059   0.0060    0.0087    0.0187



    
    
