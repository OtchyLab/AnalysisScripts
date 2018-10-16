%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes spectrogram and syllable duration histograms
%for selected days of the dataset.

%% Set function constants

%File location
mother = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis';
file = 'LW60_0606-0628_intervals.mat'; %<=== Update this to load different file

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
cmap = load('tCAF_cmap.mat');

%Load the intervals file (from Metermeter2 output)
fLoc = [mother, filesep, file];
load(fLoc, 'pData');

%% Define the figure to make
figure(99); clf
set(gcf, 'Units', 'inches', 'Position', [10.25, 4, 9.75, 9.25])


%% Plot spectrograms of the first and last days in the dataset

%Pointer to which days and files to plot ([first day, last day])
% days = [1, 2];% Baseline example
% days = [6, 12];% Drive up example
days = [12, 17];% Drive down example
pntr = [1, 1];

%Plot spectrograms
r = 1;
for i = days %days
    %Prep signal
    sig = pData(i).audio{pntr(r)};
    filtSig = filtfilt(BP_b, BP_a, sig);
    
    subplot(4, 2, r)
    displaySpecgramQuick(filtSig, fs,[300,8000],[],0);
    h(r) = gca;
    
    %Format axis
    set(h(r), 'Box', 'off', 'TickDir', 'out')
    title(['Motif on Day ' num2str(i)])
    r = r+1;
end
h(1).Colormap = cmap.cmap;
h(2).Colormap = cmap.cmap;
%% Plot the histograms of the duration of each syllable

%Combine gap intervals (if requested)
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

%Calculate the histograms for each syllable each day


%Plot histograms (D1)
for i = 1:numInts
    subplot(4, 2, (2*i+1)); cla
    h = histogram(squeeze(data(1,:,i)), 'LineWidth', 1, 'FaceColor', 'b', 'FaceAlpha', 0.6);
    hold on
    
    %Format
    h.BinLimits = rg(i,:);
    h.BinWidth = binSize;
    set(gca, 'Box', 'off', 'TickDir', 'out') 
    ylim([0, 35])
    ylabel(['Syllable ' num2str(i) ' Counts'])
end
xlabel('Interval Duration (ms)')

%Plot histograms (D2)
for i = 1:numInts
    subplot(4, 2, (2*i+2))
    h = histogram(squeeze(data(2,:,i)), 'LineWidth', 1, 'FaceColor', 'r', 'FaceAlpha', 0.6);

    %Format
    h.BinLimits = rg(i,:);
    h.BinWidth = binSize;
    set(gca, 'Box', 'off', 'TickDir', 'out')
    ylim([0, 35])
end
xlabel('Interval Duration (ms)')







