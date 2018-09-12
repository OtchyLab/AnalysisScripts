function S = featureBuilder(snip)
%This script is for analyzing data from the fictive singing experiments.
%The purpose of this script is to build a feature set from each audio snip
%as extracted using the makeStimSnips.m script.
%
%INPUTS
%   snip = 1xN array containing the audio timeseries for single trial
%          (sampled at 44150Hz; initial filtering done prior)
%
%OUTPUTS
%   featureMatrix = SxT matrix containing S features calculated over T 
%                   1ms timebins. the feature sets (S) are:
%                   S(1,:) = Amplitude envelop (log)
%                   S(2,:) = Frequency Mod
%                   S(3,:) = Amplitude Mod
%                   S(4,:) = Weiner entropy
%                   S(5,:) = Pitch
% Written by TMO; last mod 06/21/2018

bDebug = true;
if bDebug
    
end

%Set constants for the file
fs = 44150; %sampling rate (in Hz)

%Smoothed audio envelop
% envlp = smooth((x-mean(x)).^2, 221, 'moving')';
envlp = envelope(snip, 221, 'rms');
% [envlp, ~] = aSAP_getLogPower(snip, fs);
    
%Calculate the standard SAP timeseries
features = koenigSpectral(snip, fs);

%Strip to vars
FM = features.FM';
AM = features.AM';
entropy = features.Entropy';
good = features.PitchGoodness';


%Define the yin parameters ({} ==> will be default/calculated)
P = [];
P.minf0 = 350;  %    Hz - minimum expected F0 (default: 30 Hz)
P.maxf0 = 1500; %    Hz - maximum expected F0 (default: SR/(4*dsratio))
P.thresh = 0.1; %   threshold (default: 0.1)
P.relfag = 1;   %   if ~0, thresh is relative to min of difference function (default: 1)
P.hop = 44;      %   samples - interval between estimates (default: 32)
%P.range = [];    %   samples - range of samples ([start stop]) to process
%P.bufsize = [];  %   samples - size of computation buffer (default: 10000)
P.sr = fs;		%   Hz - sampling rate (usually taken from file header)
%P.wsize = [];	%   samples - integration window size (defaut: SR/minf0)
P.lpf = 300;		%   Hz - intial low-pass filtering (default: SR/4)
P.shift = 0;	% 	0: shift symmetric, 1: shift right, -1: shift left (default: 0)
yinPitch = [];

%Make YIN pitch measurement
Rs = yin(snip',P);
octs = Rs.f0; %Pitch measure (in octaves)
yinPitch = 440.*(2.^octs(2:end-1)); %Convert to Hz

%Timing section (necessary?)

%Copy the measurements over to the output structure (resample to the SAP feature array length)
S(1,:) = resample(envlp, length(AM), length(envlp)); %Amplitude envelop (log)
S(2,:) = FM; %Frequency Mod
S(3,:) = AM; %Amplitude Mod
S(4,:) = entropy; %Weiner entropy
S(5,:) = resample(yinPitch, length(AM), length(yinPitch)); %Pitch
S(6,:) = good;

% Smoothing of the features output
for i = 1:size(S,1)
    S2(i,:) = smooth(S(i,:)', 15, 'moving')';
    
end
S=S2;

% %Add in the timing gaussians
% t = -50:50;
% % overlap = 0.25;
% n = normpdf(t, 0, 1);
% % numT = ceil(size(S,2)/(numel(t)*(1-overlap)));
% timeS = zeros(1, size(S,2));
% timeS(1:101) = n;
% timeS(200:300) = n;
% % 
% % for i = 1:numT
% %    starts = floor(numel(t)*(1-overlap)*(i-1)+1);
% %    ends = starts+numel(t)-1;
% %    timeS(i, starts:ends) = n;
% % end
% % timeS = timeS(:,1:size(S,2));
% 
% S = [S; timeS];









