%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes spectrogram and syllable duration histograms
%for selected days of the dataset.

%% Set function constants
clear

%File location
% mother = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis';
mother = '/Users/Tim/Dropbox';
file = 'LW60_0606-0719_intervals.mat'; %<=== Update this to load different file

%Signal processing defaults
fs = 44150; %sampling rate (in Hz)
HP = 300; %High pass cutoff (in Hz)
LP = 8000; %Low pass cutoff (in Hz)
poles = 7;
gain = 10;
win = 220;
step = 44;
nfft = 512;

%Constants for Bandpass Audio
HP_fNorm = HP/(fs/2);
LP_fNorm = LP/(fs/2);
[BP_b,BP_a] = butter(poles,[HP_fNorm LP_fNorm]);

%Set colormap
cols=get(gca,'DefaultAxesColorOrder');

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, file];
load(fLoc, 'pData');


%% Plot spectrograms of the first and last days in the dataset

%Define the figure to make
figure(98); clf
set(gcf, 'Units', 'inches', 'Position', [0.500, 7.75, 7.5, 2.75])
load('tCAF_cmap2.mat');

%Pointer to which days and files to plot)
pDay = 6;% Baseline example
pntr = 1;

%Prep signal
sig = pData(pDay).audio{pntr};
filtSig = filtfilt(BP_b, BP_a, 10.*(sig-mean(sig)));

%Plot spectrograms
displaySpecgramQuick(filtSig, fs,[300,8000],[-3,23],0);
h = gca;
h.Colormap = cmap;

%Format axis
set(h, 'Box', 'off', 'TickDir', 'out')

%% Plot the histograms for the designated days of the file

%Pointer to which days to plot (can add arbitrary many)
days = [];
days(1) = 6;% Baseline day 1
days(2) = 12;% Peak drive 1
days(3) = 17;% Return to baseline 1
numDays = numel(days);

%Combine gap intervals
numInts = 3;
ints = [1,2; 3,4; 5,nan];
data = [];
r = 1;
for i = days %days
    for j = 1:numInts %intervals
        if ~isnan(ints(j,2))
            data(r,:,j) = pData(i).intervals(:,ints(j,1)) + pData(i).intervals(:,ints(j,2));
        else
            data(r,:,j) = pData(i).intervals(:, ints(j,1));
        end
        
    end
    r = r+1;
end

%Estimate the range of the data
binSize = 2; %bin size in ms
rg = [];
for i = 1:numInts
    t = data(:,:,i); %all data for each syllable
    rg(i,1) = min(t(:))-binSize;
    rg(i,2) = max(t(:))+binSize;
end

%% Calculate and plot the histograms for each syllable each day (all on one axis)

%Define the figure to make
figure(99); clf
set(gcf, 'Units', 'inches', 'Position', [0.25, 1.75, 4, 9])

%Plot histograms
h = [];
yl = [20, 30, 40];
xl = [220,300; 195,275; 140,220];
for i = 1:numInts
    q = subplot(3, 1, i); cla
    hold on
    %Plot the histo
    q.ColorOrderIndex=1;
    for j = 1:numDays
        h = histogram(squeeze(data(j,:,i)), 'LineWidth', 2, 'EdgeColor', 'auto', 'EdgeAlpha', 0.65);
        h.BinLimits = rg(i,:);
        h.BinWidth = binSize;
        h.DisplayStyle = 'stairs';
        h.Normalization = 'count';
    end
    
    q.ColorOrderIndex=1;
    for j = 1:numDays
        %Plot the mean
        m = mean(squeeze(data(j,:,i)));
        n = scatter(m, 5, 'o');
        n.MarkerEdgeColor = 'none';
        n.MarkerFaceColor = cols(j,:);
        n.SizeData = 160;
        n.MarkerFaceAlpha = 0.65;
    end
    
    %Format
    set(gca, 'Box', 'off', 'TickDir', 'out')
    set(gca, 'YTick', [0,yl(i)/2,yl(i)], 'XTick', [xl(i,1):20:xl(i,2)])
    ylim([0, yl(i)])
    xlim(xl(i,:))
    ylabel(['Syllable ' num2str(i) ' Counts'])
end
xlabel('Interval Duration (ms)')
legend('Day1', 'Day2', 'Day3')






