%New attempt at plotting the example of chronic nerve recording with the
%nanoclip. 
%
% Created by TMO 06/24/2019

clear all

%Files to stitch
mother = '/Users/Tim/Dropbox/Test Recordings';
f = 'LR81RY177_3586769489_2017_08_28_08_51_29.612.dat';
file = [mother, filesep, f];

%Load from file
[chans,fs] = getChannels(file);

%% Process and plot audio

%Copy out
raw_audio = chans(1,:);

%Butterworth, zero-phase filtering
gain = 10;
hp = 300; %HighPass (Hz)
lp = 8000; %Lowpass (Hz)

HP_fNorm = hp/(fs/2);
LP_fNorm = lp/(fs/2); 
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

% Apply filter
filt_audio = filtfilt(BP_b,BP_a,raw_audio.*gain);

figure(1); clf
[s,f,t,ps] = spectrogram(filt_audio, 220, 176, 512, fs);
spec = 10*log10(abs(ps(4:94,:)));
imagesc(spec, [-65, -10]); axis xy; colormap(jet)

%% Process the neural signals
%Copy out
raw_chans = chans(2:5,:);
working_ch = 4;

%Common-mode subtraction
cmMask = true(4,1);
cmMask(working_ch) = false;
cm = mean(raw_chans(cmMask, :),1);
cms = raw_chans(working_ch) - cm;

%Thresh
hiBool = cms > 0.2;
lowBool = cms < -0.2;
cms(hiBool) = cms(hiBool)-0.35;
cms(lowBool) = cms(lowBool)+0.35;

%Butterworth, zero-phase filtering
gain = 1;
hp = 300; %HighPass (Hz)
lp = 15000; %Lowpass (Hz)

HP_fNorm = hp/(fs/2);
LP_fNorm = lp/(fs/2); 
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

% Apply filter
filt_chan = filtfilt(BP_b,BP_a,cms.*gain);

% Calculate smoothed envelop
raw_env = filt_chan.^2;
hiBool = raw_env>0.025;
raw_env(hiBool) = raw_env(hiBool)-0.025;
%env = raw_env;
env = smooth(raw_env, 1101, 'moving');

figure (2); clf
subplot(4, 1, 1)
plot(raw_chans(working_ch, :));
ylim([-0.65, 0.65])

subplot(4, 1, 2)
plot(cms);
ylim([-0.65, 0.65])

subplot(4, 1, 3)
plot(filt_chan);
ylim([-0.65, 0.65])

subplot(4, 1, 4)
plot(env);

%% Plot the real figure
figure(4); clf

%Spectrogram
subplot(4,1,1)
imagesc(spec, [-65, -10]); 
axis tight; axis xy; colormap(jet)

set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

%Raw neural data
subplot(4,1,2)
t = (1:numel(filt_chan))./fs;
plot(t, raw_chans(working_ch, :), 'k', 'LineWidth', .5);
axis tight; ylim([-0.6, 0.6])

set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', -0.5:0.5:0.5)

%CMS & filtered neural data
subplot(4,1,3)
plot(t, filt_chan, 'k', 'LineWidth', .5);
axis tight; ylim([-0.35, 0.35])

set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', -0.3:0.3:0.3)

%Envelop
subplot(4,1,4)
plot(t, env, 'k', 'LineWidth', .5);
axis tight; ylim([0, 0.006])

set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:10, 'YTick', [0, 0.005])

set(gcf, 'Units', 'Inches', 'Position', [3, 4.25, 13.5, 6.75])






