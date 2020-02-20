%Plotting eight examples of chronic nerve recording with the nanoclips
%
% Created by TMO 06/24/2019

clear all

%Files to stitch
mother = '/Users/Tim/Dropbox/Test Recordings';
f = 'LR81RY177_170828_dataset345.mat';
file = [mother, filesep, f];

%Load from file
load(file, 'neuro');
raw = neuro.raw;
aligned = neuro.aligned;
ref = neuro.ref;
aligned_TS = neuro.aligned_TS;

clear('neuro')

%% Select most correlated traces

%Cross-correlation matrix for the aligned power
pow = aligned;
corr_pow = corrcoef(pow');
% figure(1); clf
% imagesc(corr_pow);

%Find high corr sets
num_rends = size(pow,1);
[B1,Is] = sort(corr_pow, 2);
[B2,Im] = sort(B1(:,num_rends-7));
r = Im(end-3);
idx = Is(r,num_rends-7:num_rends-1);
idx = [r, idx];

%% process and plot examples

%Butterworth, zero-phase filtering
fs = 44150;
gain = 10;
hp = 240; %HighPass (Hz)
lp = 6000; %Lowpass (Hz)

HP_fNorm = hp/(fs/2);
LP_fNorm = lp/(fs/2); 
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

figure(2); clf
for i = 1:numel(idx)
    %Process the raw data
    vt = raw{idx(i)} - ref{idx(i)};
    vt_filt = filtfilt(BP_b, BP_a, vt.*gain);
    t = 1000*(1:numel(vt_filt))/fs;
    
    subplot(numel(idx),1,i)
    plot(t, vt_filt, 'k', 'LineWidth', 1)

    axis tight; ylim([-0.55, 0.55])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', -0.5:0.5:0.5)

end
set(gca, 'XTick', [0, 500])

set(gcf, 'Units', 'Inches', 'Position', [7.25, 1.75, 6.75, 9.25])






