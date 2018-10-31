%pre-post ChABC injections analysis script for Darkwa, Nerurkar, Semu, and Otchy (2018).
%This makes a few timeseries plots for example animals.

%% Set function constants
clear

%File location
mother = 'C:\Users\Tim\Desktop\';
% mother = '/Users/Tim/Dropbox';

%ChABC into RA
% file = 'LLY77_0219-0423_intervals.mat'; %<=== Update this to load different file
% file = 'LLY72_0108-0210_intervals.mat';
file = 'LR71_0108-0210_intervals.mat';

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, file];
load(fLoc, 'pData');

%Set plotting defaults
load('tCAF_cmap2.mat');

%Timepoints
pre = [-3, -2, -1];
post = [5:10, 15, 20, 25, 30];
xs = [pre, post];

%Number of pre-injection days
numPre = numel(pre);
numPost = numel(post);

%% Prepare the data and calculate summary stats

%Number of days
numDays = numel(pData);
if numel(xs) ~= numDays
    display('Check the data... number of timepoints isn''t right?')
    return
end

%Combine intervals
if size(pData(1).intervals, 2) == 5
    % 3 syllables
    numInts = 3;
    ints = [1,2; 3,4; 5,nan];
elseif size(pData(1).intervals, 2) == 7
    % 4 syllable
    numInts = 4;
    ints = [1,2; 3,4; 5,6; 7,nan];
else
    display('Wrong interval numbers... double check it')
    return
end

%Combine gap intervals (if requested)
data = nan(numDays,100,numInts);
date = []; numTrials = [];
for i = 1:numDays %days
    numTrials(i) = size(pData(i).intervals,1);
    
    %Combine intervals
    for j = 1:numInts %intervals
        if ~isnan(ints(j,2))
            data(i,1:numTrials(i),j) = pData(i).intervals(:,ints(j,1)) + pData(i).intervals(:,ints(j,2));
        else
            data(i,1:numTrials(i),j) = pData(i).intervals(:, ints(j,1));
        end
    end
    
    %Copyout datenum and numInts
    date(i) = pData(i).date(1);
    
end

%Readable dates
readableDate = datestr(date);
ts = date-date(1);

%Calculate descriptive stats over intervals and days
means = [];
stds = [];
for i = 1:numDays %days
    %
    for j = 1:numInts %intervals
        means(i,j) = nanmean(squeeze(data(i,:,j)));
        stds(i,j) = nanstd(squeeze(data(i,:,j)),1);
    end
end
cvs = stds./means;


%% Breakout into useful units

%Mean over pre intervals
preInts = mean(means(1:numPre,:),1);
preCV = 100*mean(cvs(1:numPre,:),1);

%Extract post intervals
postMeans = means(numPre+1:end,:);
postStds = stds(numPre+1:end,:);
postCV = 100*cvs(numPre+1:end,:);

%Normalize as required
normPost = postMeans./preInts;
normPostStd = postStds./preInts;
normPostSem = normPostStd./sqrt(numTrials(numPre+1:end))';

%% Plot the output

%Define the figure to make
figure(69); clf
set(gcf, 'Units', 'inches', 'Position', [11.5, 7, 7.25, 3.25])

%Normalized interval durations
subplot(1,2,1); cla
for i = 1:numInts
    %Plot each trace
    e = plot([-1, post], [1, smooth(normPost(:,i)', 3)']); hold on
    %e.Color = [0.5, 0.5, 0.5];
    e.LineWidth = 0.5;
    e.LineStyle = '-';

end
shadedErrorBar(post, smooth(mean(normPost, 2)', 3)', std(normPost, 1, 2)'./sqrt(numInts), [], 1)

%Formatting
xlim([-5, 45]); ylim([0.8, 1.2])
xlabel('Days Post-Lesion')
ylabel('Change in Interval Duration (% Change of Pre-Lesion)')

set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', 0.8:0.2:1.2, 'YTickLabel', {-20, 0, 20})

%Normalized interval CV
subplot(1,2,2); cla
for i = 1:numInts
    %Plot each trace
    e = plot([-1, post], [preCV(:,i), smooth(postCV(:,i)', 3)']); hold on
    %e.Color = [0.5, 0.5, 0.5];
    e.LineWidth = 0.5;
    e.LineStyle = '-';

end
errorbar(-1, mean(preCV), std(preCV,1)./sqrt(numInts), 'ok')
shadedErrorBar(post, smooth(mean(postCV, 2)', 3)', std(postCV, 1, 2)'./sqrt(numInts), [], 1)

%Formatting
xlim([-5, 45]); ylim([0, 6])
xlabel('Days Post-Lesion')
ylabel('Interval Duration CV (%)')

set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', 0:2.5:6)






