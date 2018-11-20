%pre-post ChABC injections analysis script for Darkwa, Nerurkar, Semu, and Otchy (2018).
%This makes a few bar plots for example animals. The Final few sections
%allow collating datat across animals for summary plots

%% Set function constants
clear

%File location
mo = '/Users/Tim/Dropbox/Summary Ints';
outFile = 'longTemporalData_v1.mat';

%ChABC in HVC
type = 'HVC_ChABC';
mother = '/Users/Tim/Dropbox/Summary Ints/HVC';
% file = 'LLY72_0108-0210_intervals.mat';
% file = 'LR71_0108-0210_intervals.mat';
% file = 'LR47_0216-0321_intervals.mat';
% file = 'LR58_0216-0321_intervals.mat';
file = 'LW58_0717-0819_intervals.mat';
% file = 'LW60_0717-0819_intervals.mat';

%ChABC into RA
% type = 'RA_ChABC';
% mother = '/Users/Tim/Dropbox/Summary Ints/RA';
% file = 'LLY77_0219-0423_intervals.mat'; %<=== Update this to load different file
% file = 'LR60_0109-0212_intervals.mat';
% file = 'LLR23_0109-0212_intervals.mat';
% file = 'LLY41RB_0219-0325_intervals.mat';

%PBS into HVC
% type = 'HVC_PBS';
% mother = '/Users/Tim/Dropbox/Summary Ints/PBS';
% file = 'LR41_0330-0502_intervals.mat';
% file = 'LLB88_0303-0406_intervals.mat';

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, file];
load(fLoc, 'pData');

sp = regexp(file, '_', 'split');

%Timepoints
pre = -3:-1;
p5 = 5:6;
p10 = 9:12;
p20 = 18:22;
p30 = 29:34;

%% Prepare the data and calculate summary stats

%Number of days
numDays = numel(pData);

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
ts = date-(date(3)+1);

%% Pare down the total dataset to the timepoints we care about

%Refresh var
procData = [];
j = 0;
idx = 1:size(data,1);

%Pre assembly
mask = ismember(ts, pre);
tmp = data(mask,:,:);
q = [];
for i = idx(mask)
    q = [q; squeeze(data(i,:,:))];
end
j = j+1;
procData(j).bird = sp{1};
procData(j).type = type;
procData(j).pnt = 'Pre';
procData(j).intDur = q;

%P5
mask = ismember(ts, p5);
tmp = data(mask,:,:);
q = [];
for i = idx(mask)
    q = [q; squeeze(data(i,:,:))];
end
j = j+1;
procData(j).pnt = 'P5';
procData(j).intDur = q;

%P10
mask = ismember(ts, p10);
tmp = data(mask,:,:);
q = [];
for i = idx(mask)
    q = [q; squeeze(data(i,:,:))];
end
j = j+1; 
procData(j).pnt = 'P10';
procData(j).intDur = q;

%P20
mask = ismember(ts, p20);
tmp = data(mask,:,:);
q = [];
for i = idx(mask)
    q = [q; squeeze(data(i,:,:))];
end
j = j+1;
procData(j).pnt = 'P20';
procData(j).intDur = q;

%P30
mask = ismember(ts, p30);
tmp = data(mask,:,:);
q = [];
for i = idx(mask)
    q = [q; squeeze(data(i,:,:))];
end
j = j+1;
procData(j).pnt = 'P30';
procData(j).intDur = q;

%% Process each timepoint into useful quantities

%Number of timepoints
numBins = numel(procData);

%Sequentially process data
for i = 1:numBins
    %Motif durations
    procData(i).motifDur = nansum(procData(i).intDur,2);
    procData(i).motifDur(procData(i).motifDur==0) = [];
    
    %Motif descriptive stats
    procData(i).motifMean = nanmean(procData(i).motifDur);
    procData(i).motifStd = nanstd(procData(i).motifDur);
    procData(i).motifSem = nanstd(procData(i).motifDur)/sqrt(numel(procData(i).motifDur));

    %Interval descriptive stats
    procData(i).intMean = nanmean(procData(i).intDur,1);
    procData(i).intStd = nanstd(procData(i).intDur,1);
    procData(i).intStd = nanstd(procData(i).intDur,1)/sqrt(numel(procData(i).motifDur));
    procData(i).intCV = 100.*nanstd(procData(i).intDur,1)./nanmean(procData(i).intDur,1);
    procData(i).intMCV = mean(procData(i).intCV);
    procData(i).intSCV = std(procData(i).intCV, 1);
    
    if isempty(procData(i).intCV)
        procData(i).intCV = NaN;
    end
end


%% Plot the output

%Define the figure to make
h = figure(70); clf
set(gcf, 'Units', 'inches', 'Position', [5.25,6.75,9,3.5])

%Arrange for plotting
mMeans = getFieldVector(procData, 'motifMean');
mStd = getFieldVector(procData, 'motifStd');
iCV = cell2mat(getFieldVectorCell(procData, 'intCV')');
iMCV = getFieldVector(procData, 'intMCV');
iSCV = getFieldVector(procData, 'intSCV');
mLabels = getFieldVectorCell(procData, 'pnt');
xs = 1:numel(mMeans);

%Motif duration mean and variance
subplot(1,2,1); cla

b = bar(xs, mMeans); hold on
b.LineWidth = 1.5;
b.FaceColor = [0.5, 0.5, 0.5];

e = errorbar(xs, mMeans, mStd, 'Color', 'k', 'LineStyle', 'none');
e.LineWidth = 1.5;

%Formatting
xlim([0, 6]); ylim([0, 1000])
ylabel('Motif Length (ms)')

set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', 0:500:1000, 'XTick', xs, 'XTickLabel', mLabels)


%Motif duration mean and variance
subplot(1,2,2); cla

b = bar(xs, iMCV); hold on
b.LineWidth = 1.5;
b.FaceColor = [0.5, 0.5, 0.5];

e = errorbar(xs, iMCV, iSCV, 'Color', 'k', 'LineStyle', 'none'); hold on
e.LineWidth = 1.5;

n = size(iCV, 2);
s = plot(repmat(xs,n,1)', iCV, 'o');
for i = 1:numel(s)
    s(i).MarkerSize = 10;
end

%Formatting
xlim([0, 6]); ylim([0, 12])
ylabel('Interval CV (%)')

set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', 0:5:10, 'XTick', xs, 'XTickLabel', mLabels)

savefig(h, [mother, filesep, file(1:end-13), 'tempStats.fig'])

% %% Pump it to the running save file
% 
% %%%%%%%%%%%%%%%%%%%%%
% % Save to output
% %%%%%%%%%%%%%%%%%%%%%
% outName = [mo, filesep, outFile];
% m = exist(outName);
% if m == 2 
%     %File already exists
%     load(outName, 'stats')
%     stats(end+1) = procData;
% else
%     %No file yet created
%     stats = procData;
% end
% 
% %Save the updated data to file
% save(outName, 'stats')
% 
% display('done with the proc file')

%% Strip the good stuff and save to file

summaryOut = 'summaryStats_all_v1.mat';

out = [];
out.bird = procData(1).bird;
out.type = procData(1).type;
out.bins = getFieldVectorCell(procData,'pnt');
out.motifMean = getFieldVector(procData,'motifMean');
out.CVMean = getFieldVector(procData,'intMCV');

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
outName = [mo, filesep, summaryOut];
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

display('done with the summary file')

%% Plot the summary across animals

%Add the difference matrix to the plot
h1 = figure('Name', 'Motif Length Change'); clf
set(gcf, 'Units', 'Inches', 'Position', [11.75,3.25,5.75,4.5])

%Get indices
tVect = getFieldVectorCell(stats, 'type');
typelist = unique(tVect);
typelist = typelist([1, 3, 2]);
xs = 1:5;

mmCell = getFieldVectorCell(stats, 'motifMean');
mmVect = cell2mat(mmCell');
for i = 1:numel(typelist)
    %Masking
    mask = ismember(tVect, typelist(i));
    subset = mmVect(mask,:);
    
    %Normalize by pre
    nsubset = subset./subset(:,1);
    
    %Descriptive stats
    ms = nanmean(nsubset, 1);
    ss = nanstd(nsubset, 1, 1);

    %Scatter totals
    r = repmat(xs+6*(i-1), size(nsubset,1), 1);
    xss = r(:);
    yss = nsubset(:); 
%     yss = nsubset(:);
%     xss = repmat(xs+6*(i-1), 1, size(nsubset,1));
%     
    %Plot the bar
    b = bar(xs+6*(i-1), ms, 'k'); hold on
    b.LineWidth = 1.5;
    b.FaceColor = [0.5, 0.5, 0.5];

    e = errorbar(xs+6*(i-1), ms, ss, 'Color', 'k', 'LineStyle', 'none');
    e.LineWidth = 1.5;

    s = scatter(xss, yss, 'ok', 'jitter','on', 'jitterAmount',0.05);
    s.SizeData = 100;
    
end
xlim([0, 18]); ylim([0.5, 1.5])
ylabel('Normalized Change in Motif Length')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [1:5, 7:11, 13:17], 'XTickLabel', stats(1).bins)
set(gca, 'YTick', 0:.5:1.5)

% savefig(h, [mos, filesep, subset(2).filename(1:end-15), '.fig'])



h2 = figure('Name', 'CV Change'); clf
set(gcf, 'Units', 'Inches', 'Position', [11.75,3.25,5.75,4.5])

mmCell = getFieldVectorCell(stats, 'CVMean');
mmVect = cell2mat(mmCell');
for i = 1:numel(typelist)
    %Masking
    mask = ismember(tVect, typelist(i));
    subset = mmVect(mask,:);
    
    %Descriptive stats
    ms = nanmean(subset, 1);
    ss = nanstd(subset, 1, 1);

    %Scatter totals
    r = repmat(xs+6*(i-1), size(subset,1), 1);
    xss = r(:);
    yss = subset(:);
    
    %Plot the bar
    b = bar(xs+6*(i-1), ms, 'k'); hold on
    b.LineWidth = 1.5;
    b.FaceColor = [0.5, 0.5, 0.5];

    e = errorbar(xs+6*(i-1), ms, ss, 'Color', 'k', 'LineStyle', 'none');
    e.LineWidth = 1.5;

    s = scatter(xss, yss, 'ok', 'jitter','on', 'jitterAmount',0.05);
    s.SizeData = 100;
    
end
xlim([0, 18]); ylim([0, 12])
ylabel('Song Interval CV (%)')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [1:5, 7:11, 13:17], 'XTickLabel', stats(1).bins)
set(gca, 'YTick', [0,12])


