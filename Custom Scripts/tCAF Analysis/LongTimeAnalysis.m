%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes timeseries plots of the syllable durations
%and changes. Also saves some summary data to a file to group data over
%birds and conditions.

%% Set function constants
clear 

%File location
mother = 'C:\Users\Tim\Desktop';

%Pointer to which days and files to extract (LW60 Control)
out.file = 'LW60_0606-0719_intervals.mat'; %<=== Update this to load different file
out.batch = 'Control';
days(1) = 6;% Baseline day 1
days(2) = 12;% Peak drive 1
days(3) = 17;% Return to baseline 1
days(4) = 24;% Baseline day 2
days(5) = 28;% Peak drive 2
days(6) = 40;% Spontaneous return day

% %Pointer to which days and files to extract (LW60 ChABC)
% out.file = 'LW60_0801-0919_intervals.mat'; %<=== Update this to load different file
% out.batch = 'ChABC';
% days(1) = 6;% Baseline day 1
% days(2) = 12;% Peak drive 1
% days(3) = 17;% Return to baseline 1
% days(4) = 24;% Baseline day 2
% days(5) = 28;% Peak drive 2
% days(6) = 36;% Spontaneous return day

% %Pointer to which days and files to extract (LW58 Control)
% out.file = 'LW58_0606-0719_intervals.mat'; %<=== Update this to load different file
% out.batch = 'Control';
% days(1) = 6;% Baseline day 1
% days(2) = 13;% Peak drive 1
% days(3) = 17;% Return to baseline 1
% days(4) = 22;% Baseline day 2
% days(5) = 27;% Peak drive 2
% days(6) = 29;% Spontaneous return day

% %Pointer to which days and files to extract (LW58 ChABC)
% out.file = 'LW58_0801-0919_intervals.mat'; %<=== Update this to load different file
% out.batch = 'ChABC';
% days(1) = 3;% Baseline day 1
% days(2) = 4;% Peak drive 1
% days(3) = 7;% Return to baseline 1
% days(4) = 9;% Baseline day 2
% days(5) = 10;% Peak drive 2
% days(6) = 12;% Spontaneous return day

% %Pointer to which days and files to extract (LY80 Control)
% out.file = 'LY80_0606-0719_intervals.mat'; %<=== Update this to load different file
% out.batch = 'Control';
% days(1) = 2;% Baseline day 1
% days(2) = 6;% Peak drive 1
% days(3) = 12;% Return to baseline 1
% days(4) = 22;% Baseline day 2
% days(5) = 26;% Peak drive 2
% days(6) = 34;% Spontaneous return day

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, out.file];
load(fLoc, 'pData');

%% Plot long time series of syllable durations

%Bird specific details
if size(pData(1).intervals, 2) == 5
    % %LW60 & LY80
    numInts = 3;
    ints = [1,2; 3,4; 5,nan];
elseif size(pData(1).intervals, 2) == 7
    % %LW58
    numInts = 4;
    ints = [1,2; 3,4; 5,6; 7,nan];
else
    return
end

%How many days in this dataset
numDays = numel(pData);

%Combine gap intervals (if requested)
data = nan(numDays,100,numInts);
date = [];
for i = 1:numDays %days
    numTrials = size(pData(i).intervals,1);
    
    %Combine intervals
    for j = 1:numInts %intervals
        if ~isnan(ints(j,2))
            data(i,1:numTrials,j) = pData(i).intervals(:,ints(j,1)) + pData(i).intervals(:,ints(j,2));
        else
            data(i,1:numTrials,j) = pData(i).intervals(:, ints(j,1));
        end
    end
    
    %Copyout datenum
    date(i) = pData(i).date(5);
end

%Readable dates
readableDate = datestr(date);
ts = date-date(1);
%Calculate descriptive stats
means = [];
stds = [];
for i = 1:numDays %days
    %Combine intervals
    for j = 1:numInts %intervals
        means(i,j) = nanmean(squeeze(data(i,:,j)));
        stds(i,j) = nanstd(squeeze(data(i,:,j)),1);
    end
end


%% Setup figure
figure(100); clf
set(gcf, 'Units', 'inches', 'Position', [10.25, 3.5, 6.5, 9.5])

%Plot the syllable durations
subplot(2,1,1)
for i = 1:numInts
    errorbar(ts, means(:,i), stds(:,i), 'Marker', '.', 'MarkerSize', 10, 'LineWidth', 1.5); hold on
end

%Format
set(gca, 'Box', 'off', 'TickDir', 'out')
ylabel('Syllable Duration (ms)')
xlabel('Time (days)')

%Plot the change in syllable duration
baseDays = 1:7;
meanBase = mean(means(baseDays,:), 1);

subplot(2,1,2)
diffs = means-meanBase;
col = {'b', 'r', 'g'};
for i = 1:numInts
%     errorbar(ts, diffs(:,i), stds(:,i), 'Marker', '.', 'MarkerSize', 10, 'LineWidth', 1.5); hold on
    shadedErrorBar(ts, diffs(:,i)', stds(:,i)'./10, col{i}, 1); hold on
end
xl = xlim;
line(xl, [0,0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5)

%Format
set(gca, 'Box', 'off', 'TickDir', 'out')
ylabel('Change in Duration (ms)')
xlabel('Time (days)')

%% Spotlight plots
figure(101); clf
set(gcf, 'Units', 'inches', 'Position', [10.25, 3.5, 6.5, 4.75])
col = {'b', 'r', 'g'};

%Plot the Active Drive
subplot(1,2,1)
baseDays1 = 5:6;
driveDays1 = 5:17;
meanBase1 = mean(means(baseDays1,:), 1);
stds1 = stds(driveDays1,:);
diffs1 = means(driveDays1,:)-meanBase1;
for i = 1:numInts
    shadedErrorBar((1:numel(driveDays1))-numel(baseDays1), diffs1(:,i)', stds1(:,i)'./10, col{i}, 1); hold on
end
xlim([-3, 12])
ylim([-5, 25])
xl = xlim;
line(xl, [0,0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5)

%Format
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', -5:10:25)
ylabel('Change in Duration (ms)')
xlabel('Time (days)')

%Plot the Spontaneous Drive
subplot(1,2,2)
baseDays2 = 22:23;
driveDays2 = 22:43;

meanBase2 = mean(means(baseDays2,:), 1);
stds2 = stds(driveDays2,:);
diffs2 = means(driveDays2,:)-meanBase2;
for i = 1:numInts
    shadedErrorBar((1:numel(driveDays2))-numel(baseDays2), diffs2(:,i)', stds2(:,i)'./10, col{i}, 1); hold on
end
xlim([-3, 21])
ylim([-5, 25])
xl = xlim;
line(xl, [0,0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5)

%Format
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', -5:10:25)
ylim([-5, 25])
xlabel('Time (days)')

%% Extract the interval durations at critical timepoints

%Suck out those days
out.means_reduced = means(days,:);
out.stds_reduced = stds(days,:);
out.diff_reduced = diffs(days,:);

%Get the quantities to plot
% out.deltas = diff(out.diff_reduced,1);
out.deltas = diff(out.means_reduced,1);
out.driveLength = diff(date(days));
out.shiftRate = out.deltas./out.driveLength';

%Group
m = false(size(out.shiftRate, 2),1);
m(end-1) = true;
out.tStretch = out.deltas([1,4],m);
out.ntStretch = out.deltas([1,4], ~m);
out.rateUp = out.shiftRate([1,4],m);
out.rateDown = out.shiftRate(2,m);
out.rateSpont = out.shiftRate(5,m);

%%
%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
outName = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis\summaryStats_v2.mat';
m = exist(outName);
if m == 2 
    %File already exists
    load(outName, 'stats')
    stats(end+1) = out;
else
    %No file yet created
    stats = out;
end

%Save the updated data to file
save(outName, 'stats')

display('done')















