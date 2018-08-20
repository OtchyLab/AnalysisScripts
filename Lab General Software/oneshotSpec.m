%The purpose of this script is to generate a reasonably nice looking spectrogram from a TMO-formatted *.wav file. Key
%modification points in the code are:
% filename -- the name of the file to load and display -- either absolute or relative path
% clims -- scaling values to use for the final image. You can get a sense of this from the histogram in Figure 1.

%%
%load data from file
%40dph
filename = 'C:\Users\Tim\Desktop\JGO Specs\040dph\Pur304_3362155583_2010_07_16_16_06_22.wav';

%60dph
% filename = 'C:\Users\Tim\Desktop\JGO Specs\60dph\Sil167_3344505891_2009_12_24_08_24_50.wav';

%100dph
% filename = 'C:\Users\Tim\Desktop\JGO Specs\100dph\Sil167_3348264211_2010_02_05_20_23_30.wav';

[rawaudio, ~] = audioread(filename);
fs = 44150;

%Do some basic bandpass filtering
if ~isempty(rawaudio)
    %Constants for Bandpass Audio (300-6500kHz)
    HP_fNorm = 300/(fs/2);
    LP_fNorm = 6500/(fs/2);
    [BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
    audio = filtfilt(BP_b,BP_a,rawaudio);
end

%%
% Generate basic spectrogram
[S,F,T,P] = spectrogram((audio/(sqrt(mean(audio.^2)))),220,220-44,512,fs);
strt = 4; stp = 94; %These bins correspond to 300-8000Hz
spec = abs(P(strt:stp,:));

figure(1)
hist(log10(spec(:)),100)

figure(2)
%For 40dph
clims = [-5.75, -.75];

%For 60dph
% clims = [-5.75, -.75];

%For 100dph
% clims = [-5.75, -.75];

imagesc(log10(spec), clims)
axis xy; axis tight
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', [0.5, 89.5], 'YTickLabel', {'0', '8000'})