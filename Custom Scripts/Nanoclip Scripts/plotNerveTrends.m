%Script plots the trends in neural dynamics over 30 days of recording
%for the nanoclip paper. Also computes the required significance tests as
%reported in the paper. The script reads in the rawdata from a previously
%constructed .mat file containing all of the data from each condition.
%
% Written by TMO 05-15-19

%Clear it all
clear all

%Load from file
load('/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/neuroDynamicsTrend.mat')
if ~exist('corrDat', 'var')
    disp('Uh-oh... missing data?')
    return
end

%% Plot the individual bird traces
figure(100); clf
%set(gcf,'Units', 'Inches', 'Position', [6, 5.25, 7.75, 6])

%Plot constants
days = 1:4;
birds = 1:3;

%Plot correlation data
subplot(1,2,1)
for i = birds
    plot(days, corrDat(i,:), '-ok', 'LineWidth', 1); hold on
end

%Format figure
axis square
xlim([0.5, 4.5]); ylim([0.5, 1.05])
xlabel('Time (days)'); ylabel('Correlation')
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0.5:0.5:1)


%Plot mean data
subplot(1,2,2)
for i = birds
    plot(days, meanDat(i,:), '-ok', 'LineWidth', 1); hold on
end

%Format figure
axis square
xlim([0.5, 4.5]); ylim([0.5, 1.05])
xlabel('Time (days)'); ylabel('Mean Envelop')
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0.5:0.5:1)


%% Compute significance tests

% %Group Indices
% ind = [1:3; 4:6; 7:9];
% tp = 1:numel(days);
% nidx = [1,1,1,2,2,2,3,3,3];

%Check that all data are from normal distributions (for t-tests)
corr_ksStats = []; mean_ksStats = [];
for i = 1:4
        corr_ksStats(i) = kstest(corrDat(:,i));
        mean_ksStats(i) = kstest(meanDat(:,i));
end
corr_allNormal = all(corr_ksStats(:)); % corr_allNormal = true
mean_allNormal = all(mean_ksStats(:)); % mean_allNormal = true

%Perform paired t-tests of day 1 to each subsequent day of that group
corr_hSignif = []; corr_pSignif = [];
mean_hSignif = []; mean_pSignif = [];
corr_base = corrDat(:,1);
mean_base = meanDat(:,1);
for i = 1:4
        [corr_hSignif(i), corr_pSignif(i)] = ttest(corr_base, corrDat(:,i));
        [mean_hSignif(i), mean_pSignif(i)] = ttest(mean_base, meanDat(:,i));
end

a = 1;

%Hardcoded summary:
%
% corr_allNormal = true
% mean_allNormal = true
% 
% corr_hSignif =      NaN     0     0     0
% 
% corr_pSignif =      NaN    0.2948    0.0544    0.1990
%     
% mean_hSignif =      NaN     0     0     0
%     
% 
% mean_pSignif =      NaN    0.0653    0.0941    0.1110



    
    
