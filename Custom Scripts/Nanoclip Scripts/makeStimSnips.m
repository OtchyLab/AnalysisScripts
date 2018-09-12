% function [snips, starts, files] = makeStimSnips
%This script is for analyzing data from the fictive singing experiments.
%The purpose of this script is to scan the passed stm files, and parse the
%long files to extract the audio around the snip for analysis in other
%scripts.
%
%INPUTS
%   wavFiles = string or cell array of strings listing containing the
%              absolute paths of the stm files to be processed. If empty,
%              invokes propt for user selection GUI
%
%OUTPUTS
%   snips = nxm matrix of n audio snippets of m samples long
%   starts = nx1 array with an index of the sample point at which the snips
%            were extracted
%   files = nx1 cell array listing the filenames of the raw recording files
%           for each snip
%
% Written by TMO; last mod 06/21/2018
wavFiles = [];
bDebug = false;
if bDebug
    wavFiles = {'C:\Users\Tim\Desktop\Anestetized Singing Bird\2018-05-08\LW28_3608654213_2018_05_08_15_56_53.082.stm', 'C:\Users\Tim\Desktop\Anestetized Singing Bird\2018-05-08\LW28_3608654351_2018_05_08_15_59_10.741.stm'};
end

%Set constants for the file
fs = 44150; %sampling rate (in Hz)
HP = 1000; %High pass cutoff (in Hz)
LP = 8000; %Low pass cutoff (in Hz)
poles = 7;
gain = 100;
thresh = 2.5; %5VDC signal, so something lower should work fine, I think...
vthresh = 0.1;
preT = 100; %amount of sound before the stim marker to snip (in ms)
postT = 400; %amount of sound after the stim marker to snip (in ms)

%Constants for Bandpass Audio
HP_fNorm = HP/(fs/2);
LP_fNorm = LP/(fs/2);
[BP_b,BP_a] = butter(poles,[HP_fNorm LP_fNorm]);

% fo = 510; q = 50; bw = (fo/(fs/2))/q;
% [BP_b,BP_a] = iircomb(round(fs/fo),bw,'notch'); % Note type flag 'notch'
% fvtool(b,a);

%Convert buffer times to sample counts
preS = round(fs * (preT/1000));
postS = round(fs * (postT/1000));

%If the passed filelist is empty, request user input
if isempty(wavFiles)
    wavFiles = uipickfiles('FilterSpec', '*.stm', 'output', 'cell');
end


%Snip the sections out
snips = []; starts = []; files = []; numSpikes = [];
for j = 1:numel(wavFiles)
    %Parse the filename
    sp = regexp(wavFiles(j), filesep, 'split');
    
    %Extract channels from file
    [channels, ~] = getChannels(wavFiles{j});
    
    %Shorten
    x = channels(1,:);
    audio = gain * filtfilt(BP_b, BP_a, (x-mean(x)));
%     audio = gain * (x-mean(x)); %no filtering
    currents = channels(3,:);
    triggers = channels(4,:);
    
    %Threshold crossings (only look at rising edge)
    p = 2:numel(triggers);
    ups = find(triggers(p)>thresh & triggers(p-1)<thresh);
    
    %Snip out the audio
    for i = 1:numel(ups)
        start = ups(i)-preS;
        ends = ups(i) + postS;
        
        %Check in you're going to outrun the bounds of the recording
        if start<1 || ends>numel(triggers)
            
        else
            %Cut out the audio and store it in the stack
            snips = [snips; audio(start:ends)];
            starts = [starts; start];
            files = [files; sp{1}(end)];
            
            %Cut out the corresponding voltage
            vSpikes = currents(start:ends);
            
            %Count how many pulses appear in the voltage snippet
            p = 2:numel(vSpikes);
            ons = find(vSpikes(p)>vthresh & vSpikes(p-1)<vthresh);
            numSpikes(size(snips,1)) = numel(ons);
        end
        
    end
    
end

%Artifact removal
% snips = removeArtifacts(snips);

%Limit the examples/snips to those that were produced with the same number
%of stim pulses
% r = 100;
% idx = (numSpikes > (mode(numSpikes)-r)) & (numSpikes < (mode(numSpikes)+r));
% snips = snips(idx,:);
% starts = starts(idx,:);
% files = files(idx,:);

%Plot for inspection
specStack = [];
filt = 4:94;
tx = [0.5:200:1000.5];
xPlots = 4;
yPlots = 4;
curPlot = 1;
curFig = 1;
startFig = 49;
for i= 1:numel(starts)
    if curPlot == 1
       h(curFig)=figure(startFig + curFig); clf
       curFig = curFig+1;
    end
    
    %Select subplot
    subplot(yPlots,xPlots,curPlot)
    
    %Make and plot "clean" spectrogram
    x = snips(i,:);
    %displaySpecgramQuick(x, fs, [0, 8000])
    [~,~,~,P] = spectrogram((x/(sqrt(mean(x.^2)))),220,220-44,512,44150);
    spec = log10(1+abs(P(filt,:)));
    imagesc(spec); colormap(jet); axis xy
    xlim([0, preT+postT+0.5]); ylim([0,90])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', tx, 'XTickLabels', floor(tx-preT), 'YTick', [0,90], 'YTickLabels', [0,8000])
    ylabel(['Plot Number ' num2str(i)])
    set(gcf, 'Units', 'inches', 'Position', [8,3,9,7])
    
    %Increment subplot indexing
    curPlot = curPlot+1;
    if curPlot > (yPlots*xPlots)
        curPlot = 1;
    end
    
end

%Request user input on which to retain
picks = input('Which syllables to capture?');
%Manual picks
% picks = [2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16, 17]; %Pattern 11
% picks = [3, 6, 7, 8, 12, 15, 16, 18, 19, 21]; %Pattern 12
% picks = [2, 3, 4, 5, 6, 7, 8, 10, 14, 15, 16, 19, 20, 22, 23]; %Pattern 15
% picks = [1, 2, 4, 5, 6, 9, 10, 12, 13, 14, 16, 20, 21]; %Pattern 16
% picks = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11]; %Pattern 17
% picks = [1, 2, 4, 5, 6, 9, 10, 11, 12, 14, 17, 18, 19, 20, 21]; %Pattern 18
% picks = [4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 28, 29, 30]; %Pattern 19
% picks = [1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 13, 14, 17, 19, 20, 21, 25, 26, 30, 31]; %Pattern 20

%Create binary index
px = false(size(idx));
px(picks) = true;

%Apply binary indexing logic
snips = snips(px,:);
starts = starts(px,:);
files = files(px,:);

% %Ask where to save
% mother = 'C:\Users\Tim\Desktop\Anestetized Singing Bird\Analyzed Data\2018-05-08 LW28\';
% fileLoc = uiputfile([mother '*.mat'], 'Where to save the snips file?');
% 
% %Save the output files
% save(fileLoc, 'snips', 'starts', 'files')