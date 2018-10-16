%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes timeseries plots of the syllable durations
%and changes. Also saves some summary data to a file to group data over
%birds and conditions.

%% Set function constants
clear 

%File location
mother = 'C:\Users\Tim\Desktop';

% %Pointer to which days and files to extract (LW60 Control)
% out.file = 'LW60_0606-0719_intervals.mat'; %<=== Update this to load different file
% out.batch = 'Control';
% days(1) = 6;% Baseline day 1
% days(2) = 12;% Peak drive 1
% days(3) = 17;% Return to baseline 1
% days(4) = 24;% Baseline day 2
% days(5) = 28;% Peak drive 2
% days(6) = 40;% Spontaneous return day

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
% days(2) = 12;% Peak drive 1
% days(3) = 17;% Return to baseline 1
% days(4) = 21;% Baseline day 2
% days(5) = 27;% Peak drive 2
% days(6) = 28;% Spontaneous return day

% %Pointer to which days and files to extract (LW58 ChABC)
% out.file = 'LW58_0801-0919_intervals.mat'; %<=== Update this to load different file
% out.batch = 'ChABC';
% days(1) = 2;% Baseline day 1
% days(2) = 4;% Peak drive 1
% days(3) = 7;% Return to baseline 1
% days(4) = 8;% Baseline day 2
% days(5) = 10;% Peak drive 2
% days(6) = 12;% Spontaneous return day

%Pointer to which days and files to extract (LY80 Control)
out.file = 'LY80_0606-0719_intervals.mat'; %<=== Update this to load different file
out.batch = 'Control';
days(1) = 6;% Baseline day 1
days(2) = 12;% Peak drive 1
days(3) = 17;% Return to baseline 1
days(4) = 20;% Baseline day 2
days(5) = 25;% Peak drive 2
days(6) = 31;% Spontaneous return day

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, out.file];
load(fLoc, 'pData');

%% Plot long time series of syllable durations

%Bird specific details
% %LW60
% numInts = 3;
% ints = [1,2; 3,4; 5,nan];

% %LW58
% numInts = 4;
% ints = [1,2; 3,4; 5,6; 7,nan];

%LY80
numInts = 3;
ints = [1,2; 3,4; 5,nan];

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


%Setup figure
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
for i = 1:numInts
    errorbar(ts, diffs(:,i), stds(:,i), 'Marker', '.', 'MarkerSize', 10, 'LineWidth', 1.5); hold on
end

%Format
set(gca, 'Box', 'off', 'TickDir', 'out')
ylabel('Change in Duration (ms)')
xlabel('Time (days)')
legend('Syllable 1', 'Syllable 2', 'Syllable 3', 'Syllable 4')

%% Extract and plot the interval durations at critical timepoints

%Suck out those days
out.means_reduced = means(days,:);
out.stds_reduced = stds(days,:);
out.diff_reduced = diffs(days,:);

%Get the quantities to plot
out.deltas = diff(out.diff_reduced,1);
out.driveLength = diff(date(days));
out.shiftRate = out.deltas./out.driveLength';

%Group
out.tStretch = out.deltas([1,4],2);
out.ntStretch = out.deltas([1,4],[1,3]);
out.rateUp = out.shiftRate([1,4],2);
out.rateDown = out.shiftRate(2,2);
out.rateSpont = out.shiftRate(5,2);

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
outName = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis\summaryStats.mat';
m = exist(outName);
if m ==2 
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















