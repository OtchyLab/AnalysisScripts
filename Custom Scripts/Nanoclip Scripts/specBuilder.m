function specOut = specBuilder(snip)
%This script is for analyzing data from the fictive singing experiments.
%The purpose of this script is to build a high resolution spectrogram from each audio snip
%as extracted using the makeStimSnips.m script.
%
%INPUTS
%   snip = 1xN array containing the audio timeseries for single trial
%          (sampled at 44150Hz; initial filtering done prior)
%
%OUTPUTS
%   specMatrix = SxT matrix containing S features calculated over T 
%                   1ms timebins. the feature sets (S) are:

% Written by TMO; last mod 06/21/2018

bDebug = false;
if bDebug
    
end

%Set constants for the file
fs = 44150; %sampling rate (in Hz)
HP = 300; %High pass cutoff (in Hz)
LP = 8000; %Low pass cutoff (in Hz)
poles = 7;
gain = 10;
win = 220;
step = 44;
nfft = 512;
preT = 95; %amount of sound before the stim marker to snip (in ms)
durT = 110; %amount of sound after the stim marker to snip (in ms)
% preT = 100; %amount of sound before the stim marker to snip (in ms)
% durT = 100; %amount of sound after the stim marker to snip (in ms)

%Constants for Bandpass Audio
HP_fNorm = HP/(fs/2);
LP_fNorm = LP/(fs/2);
[BP_b,BP_a] = butter(poles,[HP_fNorm LP_fNorm]);


%Convert buffer times to sample counts
% preS = round(fs * (preT/1000));
% postS = round(fs * (postT/1000));

%Process the audio
proAud = filtfilt(BP_b, BP_a, (snip-mean(snip)));
normAud = proAud/(sqrt(mean(proAud.^2)));

%Generate spectrograms for each rendition
strt = 4; stp = 94; %These bins correspond to ~300-8000Hz

% Calculate STFT features for both sounds (80% window overlap)
[S,f,t,P] = spectrogram(normAud,win,win-step,nfft,fs); %Alternatively, we could use the PSD...

%Trim filtered out frequency bands
SpecCube = abs(S(strt:stp,:));

%log transform
logSpec = gain*log10(SpecCube + 0.02);

%Trim in time
specOut = logSpec(:, preT:(preT+durT));

if bDebug
    figure; imagesc(t,f,specOut); axis xy
end
