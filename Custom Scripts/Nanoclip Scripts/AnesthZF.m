%This script is for analyzing data from the anesthetized singing bird
%experiment.

clear 

%Set threshold for TTL detection
fs = 44150; %sampling rate (in Hz)
thresh = 2.5; %5VDC signal, so something lower should work fine, I think...
vthresh = 0.1;
preT = 200; %amount of sound before the stim marker to snip (in ms)
postT = 400; %amount of sound after the stim marker to snip (in ms)
% focus = [3, 4, 8, 9, 11, 13]; %Pattern1
% focus = [3, 4, 6, 9, 11]; %Pattern2
% focus = [5, 7, 9, 11]; %Pattern3
% focus = [1, 2, 4, 5, 7, 9, 15, 16]; %Pattern4
% focus = [3, 6, 8, 16, 18]; %Pattern5
% focus = [1, 3, 4, 5, 6, 8, 9, 11, 12]; %Pattern6
% focus = [2, 7, 10, 12, 14, 15, 17, 22]; %Pattern7
focus = [2, 7, 8, 10, 11, 12, 17, 18, 19]; %Pattern8


%Convert buffer times to sample counts
preS = round(fs * (preT/1000));
postS = round(fs * (postT/1000));

%Constants for Bandpass Audio
HP_fNorm = 300/(fs/2);
LP_fNorm = 8000/(fs/2);
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

%Select and read in stm file
files = uipickfiles('FilterSpec', '*.stm', 'output', 'cell');

%Snip the sections out
snips = []; numSpikes = [];
for j = 1:numel(files)
[channels, ~] = getChannels(files{j});

%Shorten
x = channels(1,:);
audio = filtfilt(BP_b, BP_a, (x-mean(x)));
voltage = channels(2,:);
currents = channels(3,:);
triggers = channels(4,:);

%Threshold crossings (only look at rising edge)
p = 2:numel(triggers);
ups = find(triggers(p)>thresh & triggers(p-1)<thresh);

%Snip out the 
for i = 1:numel(ups)
    start = ups(i)-preS;
    ends = ups(i) + postS;
    
    %Check in you're going to outrun the bounds of the recording
    if start<1 || ends>numel(triggers)
        
    else
        %Cut out the audio and store it in the stack
        snips = [snips; audio(start:ends)];
        
        %Cut out the corresponding voltage
%         vSpikes = voltage(start:ends);
        vSpikes = currents(start:ends);
        
        %Count how many pulses appear in the voltage snippet
        p = 2:numel(vSpikes);
        ons = find(vSpikes(p)>vthresh & vSpikes(p-1)<vthresh);
        numSpikes(size(snips,1)) = numel(ons);
    end
    
end

end
%Limit the examples/snips to those that were produced with the same number
%of stim pulses
r = 100;
idx = (numSpikes > (mode(numSpikes)-r)) & (numSpikes < (mode(numSpikes)+r));
snips = snips(idx,:);
numSnips = size(snips,1);

%Plot for inspection
specStack = [];
figure(22); clf
filt = 4:94;
tx = [0.5:200:1000.5];
xPlots = 4;
yPlots = 4;
curPlot = 1;
curFig = 1;
startFig = 19;
for i= 1:numSnips
    if curPlot == 1
       h(curFig)=figure(startFig + curFig); clf
       curFig = curFig+1;
    end
    
    %Select subplot
    subplot(yPlots,xPlots,curPlot)
    
    %Make and plot "clean" spectrogram
    x = snips(i,:);
    %displaySpecgramQuick(x, fs, [0, 8000])
    [S,F,~,P] = spectrogram((x/(sqrt(mean(x.^2)))),220,220-44,512,44150);
    spec = log10(1+abs(P(filt,:)));
    imagesc(spec); colormap(jet); axis xy
    xlim([0, preT+postT+0.5]); ylim([0,90])
% 	ylabel('Frequency (Hz)')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', tx, 'XTickLabels', floor(tx-preT), 'YTick', [0,90], 'YTickLabels', [0,8000])
    
    if curPlot == 1
        title('Example Spectrograms')
    elseif curPlot > (yPlots-1)*xPlots && curPlot <= (yPlots*xPlots)
        xlabel('Time (ms)')
    end
    
    if curPlot == 1 || curPlot == 5 || curPlot == 9 || curPlot == 13
        ylabel('Frequency (Hz)')
    end
    
    set(gcf, 'Units', 'inches', 'Position', [8,3,14,8])
    
    %For averaging the spectrograms for plotting
    specStack(i,:,:) = spec;
    
    %Increment subplot indexing
    curPlot = curPlot+1;
    if curPlot > (yPlots*xPlots)
        curPlot = 1;
    end
    
end

%Plot the mean spectrogram
k=figure(50); clf
imagesc(squeeze(mean(specStack,1))); colormap(jet); axis xy
xlim([0, preT+postT+0.5]); ylim([0,90])
xlabel('Time (ms)'); ylabel('Frequency (Hz)')
title('Mean Spectrogram')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', tx, 'XTickLabels', floor(tx-preT), 'YTick', [0,90], 'YTickLabels', [0,8000])
set(gcf, 'Units', 'inches', 'Position', [13,8,5,3])

%Calculate and plot feature timeseries
envlp = [];
entropy = [];
pitch = [];
good = [];
am = [];
for i = 1: numSnips
    %Local copy
    x = snips(i,:);
    
    %Smoothed audio envelop
    %     envlp(i,:) = smooth((x-mean(x)).^2, 221, 'moving')';
    %     envlp(i,:) = envelope(x, 221, 'rms');
    envlp(i,:) = sum(squeeze(specStack(i,:,:)),1);
%     [envlp(i,:), ~] = aSAP_getLogPower(x, fs);
    
    %Extract spectral features
    features = koenigSpectral(x, fs);
    
    %Strip to vars
    %FM(i,:) = features.FM';
    am(i,:) = features.AM';
    entropy(i,:) = features.Entropy';
    good(i,:) = features.PitchGoodness';
    pitch(i,:) = features.Pitch_chose';
    
end
pitch(good<50) = NaN;
pitch(pitch>1000) = NaN;

%Plot feature projections
l=figure(51); clf

%Envelop
subplot(4,1,1)
plot(linspace(1,size(am,2),size(envlp,2))-preT, envlp);
axis tight;
xlim([-preT,postT])
ys = ylim; ylim([0, ys(2)*1.1]);
ylabel('Envelop (V^2)')
title('Aligned Feature Projections')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -200:200:600)

subplot(4,1,2)
plot((1:size(am,2))-preT, am);
axis tight;
xlim([-preT,postT])
% ys = ylim; ylim([0, ys(2)*1.1]);
ylabel('Amp Mod (ms^-1)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -200:200:600)

subplot(4,1,3)
plot((1:size(am,2))-preT, entropy);
axis tight;
xlim([-preT,postT])
ylim([-4, 0]);
ylabel('log Entropy')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -200:200:600)

subplot(4,1,4)
plot((1:size(am,2))-preT, pitch);
% plot(linspace(1,size(yinPitch,2),size(envlp,2))-preT, yinPitch);
axis tight;
xlim([-preT,postT])
xlabel('Time (ms)')
ylim([0, 1000]);
ylabel('Pitch (Hz)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -200:200:600)

set(gcf, 'Units', 'inches', 'Position', [3,3,5,8])

%Switch and code for saving the figures
bSaveIt = false;
if bSaveIt
    [file,path] = uiputfile('*.fig', 'Where to save the figures?');
    savefig([h, k, l], [path,file])
end

%Fine analysis of pitch for selected trials

%Define the yin parameters ({} ==> will be default/calculated)
P = [];
P.minf0 = 400;  %    Hz - minimum expected F0 (default: 30 Hz)
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

%Cycle through and extract yin-pitches
if ~isempty(focus)
    n = 1;
    for i = focus
        %Local copy
        x = snips(i,:);
        
        %Pitch measurements
        Rs = yin(x',P);
        if isempty(yinPitch)
            yinPitch = Rs;
        else
            yinPitch(n) = Rs;
        end
        %Smooth the power envelop
        xmat = windowTS(yinPitch(n).pwr, 5, 1, 'pad', 'boxcar');
        nM = nanmean(xmat,2)';
        
        %Threshold the envelop
        mEnv = nanmean(nM); sEnv = nanstd(nM,1);
        threshEnv = mEnv + sEnv;
        
        p = 2:numel(nM);
        ons = find(nM(p)>threshEnv & nM(p-1)<threshEnv);
        offs = find(nM(p)<threshEnv & nM(p-1)>threshEnv);
        
        ons = ons(ons>=200 & ons<=310); ons = ons(1);
        offs = offs(offs>=200 & offs<=310); offs = offs(end);
        %offs = min([235+n, offs]);
        octs = yinPitch(n).f0(ons:offs);
        pitches = 440.*(2.^octs);

        
        %Plot it
        figure(999)
        plot(ons:offs, pitches); hold on
%         plot(ons:offs, 440.*(2.^yinPitch(n).f0(ons:offs))); hold on
        
        %Increment pointer
        n = n+1;
        
    end

xlim([200, 310])
ylabel('Pitch (Hz)')
xlabel('Time (ms)')
set(gca, 'Box', 'off', 'TickDir', 'out')

end