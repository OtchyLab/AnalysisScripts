%tCAF analysis script for Darkwa, Nerurkar, and Otchy (2018). This
%script generates relevant stats for tCAF experiments that can be 
%collapsed across animals. Save the output in a running file for plotting 
%by XXX.m
%
%Created by TMO 10/26/18

%% Setup script constants
clear 

%File locations
mother = 'C:\Users\Tim\Desktop';
outName = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis\cleanedSummaryStatsv2.mat';

%Pointer to which days and files to extract (LW60 Control)
% out.file = 'LW60_0606-0719_intervals.mat'; %<=== Update this to load different file
% out.bird = 'LW60';
% out.batch = 'Control';
% base = [4:6; 21:23]; 
% dUp = [12; 26];
% dDown = 17:19;
% sDown = 40:43;

%Pointer to which days and files to extract (LW60 ChABC)
% out.file = 'LW60_0801-0919_intervals.mat'; %<=== Update this to load different file
% out.bird = 'LW60';
% out.batch = 'ChABC';
% base = [2:4; 20:22]; 
% dUp = [8; 27];
% dDown = 16:18;
% sDown = 37:39;

%Pointer to which days and files to extract (LW58 Control)
% out.file = 'LW58_0606-0719_intervals.mat'; %<=== Update this to load different file
% out.bird = 'LW58';
% out.batch = 'Control';
% base = [4:6; 20:22]; 
% dUp = [13; 27];
% dDown = 17:19;
% sDown = 29:31;

%Pointer to which days and files to extract (LW58 ChABC)
out.file = 'LW58_0801-0919_intervals.mat'; %<=== Update this to load different file
out.bird = 'LW58';
out.batch = 'ChABC';
base = [1:3; 9,9,9]; 
dUp = [4; 10];
dDown = 8;
sDown = 12;

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, out.file];
load(fLoc, 'pData');

%% Prep raw measurements

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

%% Extract the interval durations at critical timepoints

%Durations and times of interest
base1 = mean(means(base(1,:),:),1);
base2 = mean(means(base(2,:),:),1);
dUp1 = means(dUp(1),:);
dUp2 = means(dUp(2),:);
dDown1 = mean(means(dDown,:),1);
sDown1 = mean(means(sDown,:),1);

dUpTime1 = ts(dUp(1))-ts(base(1,end));
dUpTime2 = ts(dUp(2))-ts(base(2,end));
dDownTime1 = ts(dDown(1))-ts(dUp(1));
sDownTime1 = ts(sDown(1))-ts(dUp(2));

%Deltas
u1 = dUp1 - base1;
u2 = dUp2 - base2;
d1 = dDown1 - dUp1;
s1 = sDown1 - dUp2;

%Mask for the target
m = false(size(base1));
m(end-1) = true;

%Packout
bdur.target = [base1(m); base2(m)];
bdur.nontarget = [base1(~m); base2(~m)];

stretch.target = [u1(m); u2(m)];
stretch.nontarget = [u1(~m); u2(~m)];

rates.up = [u1(m)/dUpTime1, u2(m)/dUpTime2];
rates.down = d1(m)/dDownTime1;
rates.spon = s1(m)/sDownTime1;

%Copy out
out.bdur = bdur;
out.stretch = stretch;
out.rates = rates;

%% Save to output
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
