% Main Data Loading functions
clear all
%Load .dat file from disk
filename = 'V:\Grn046\Spon\2014-05-15\Grn046_3483029400_2014_05_15_16_10_00.258.dat';
fs = 44150;
[rawdata, ~] = getChannels(filename);

topLevel = regexp(filename, filesep, 'split');
lowLevel = regexp(topLevel{end}, '_', 'split');
filenum = str2num(lowLevel{2});

%Load (specific) annotation record
annotName = 'V:\Grn046\Grn046_140515_annotationSPONT.mat';
load(annotName,'elements')
filenumVect = getAnnotationVector(elements, 'filenum');
rec = find(filenumVect == filenum, 1, 'first');
annotRec = elements{rec};

% Process neural data
%Constants
chan = 4; %Select channel

%Filter the data to remove unnecessary artifacts
%Constants for Bandpass Audio (1-15kHz)
HP_fNorm = 1000/(fs/2);
LP_fNorm = 15000/(fs/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
data = filtfilt(BP_b,BP_a, rawdata(chan+1, :));

%Simple threshold crossing to detect spiketimes
idx = 2:length(data);
direction = 'pos'; %apply the spike threshold to the positive/negative side
min_n = 3; max_n = 10; step_n = 0.5;
steps = min_n:step_n:max_n;
d_m = mean(data);
d_std = std(data);
for i = 1:length(steps)
    if strcmp(direction, 'pos')
        thresh(i) = d_m + (steps(i)*d_std);
    elseif strcmp(direction, 'neg')
        thresh(i) = d_m - (steps(i)*d_std);
    end
    %Find crossing times
    ucrossings{i} = find(data(idx) > thresh(i) & data(idx-1) < thresh(i))./fs; %up-crossings
    dcrossings{i} = find(data(idx) < thresh(i) & data(idx-1) > thresh(i))./fs; %down-crossings
    
    %Count the crossings
    numUp(i) = length(ucrossings{i});
    numDown(i) = length(dcrossings{i});
end

%Select which threshold to use for subplots and rasters
selection = 5;
pntr = find(steps == selection, 1, 'first');

%Plot filtered spikes w/threshold and rasters
figure(1)
cla; hold on
h1 = plot((1:length(data))/fs, data);
line([0, length(data)/fs],[thresh(pntr), thresh(pntr)],'Color', 'r')
pos = [0.4, 0.6];
for i = 1:numUp(pntr)
    line([ucrossings{pntr}(i), ucrossings{pntr}(i)], pos, 'Color', 'k')
end
axis tight;
set(gca,'Box', 'off', 'TickDir', 'out')
hold off
title(['Voltage trace and detected spikes at V(thresh) = ' num2str(thresh(pntr))])

%Plot distribution of data w/ threshold
figure(2)
cla; hold on
subplot(2,1,1)
hist(data, 250);
y = ylim;
line([thresh(pntr), thresh(pntr)], [0, y(2)], 'Color', 'r')
set(gca,'Box', 'off', 'TickDir', 'out')
hold off

%Plot the spikes detected as function of threshold
figure(2)
subplot(2,1,2)
cla
plot(steps, numUp, 'or')
set(gca,'Box', 'off', 'TickDir', 'out')

% Parse out the excluded regions of the file
%Create the segments to exclude (based on annotation record)
buffer = 250/1000; %segmenting buffer in seconds
cutouts = [(annotRec.segFileStartTimes-buffer)', (annotRec.segFileEndTimes+buffer)'];

spikes = ucrossings{pntr};
parselength = (length(data)/fs) - sum(abs(diff(cutouts')));
for i = 1:size(cutouts,1)
    spikes = spikes(spikes<=cutouts(i,1) | spikes>=cutouts(i,2));
end

figure(1)
hold on
pos = [0.5, 0.7];
for i = 1:length(spikes)
    line([spikes(i), spikes(i)], pos, 'Color', 'g')
end

% Do some statistics on the remaining data

spikeRate  = length(spikes)/parselength;













