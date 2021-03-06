%Create the spectrogram snips for the nanoclip fictive singing t-sne
%figure. Script will read in data from a folder containing aligned
%spectrograms/syllables, format them for the t-sne analysis, and save as
%a single, bird-specific file.
%
%Written by TMO 05/16/19

%Clear the workspace
clear all

%Source and destination data
sLoc = '/Users/Tim/Dropbox/Clustering Syllables';
dLoc = '/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/FigSnips';

%Constants
maxReps = 20;
snipDur = 200; %in ms

%Define carrying vars
sylType = []; %
sylSpecs = [];

%Filelist from the source
files = dir([sLoc, filesep, '*.mat']);
numFiles = numel(files);

%Process files sequentially
for i = 1:numFiles
    %Load from file
    load(files(i).name, 'audio', 'data');
    
    %Set selection indices
    numRend = numel(audio.raw);
    idx = 1:min([maxReps, numRend]);
    
    %Bin range
    brng = data.templatesyllBreaks(1):data.templatesyllBreaks(1)+20;
    smp = linspace(1, numel(brng), snipDur);
    dataBlock = audio.aligned_audioCube(idx, :, brng);
    
    %Process renditions sequentially
    tmp = [];
    for j = idx
        for k = 1:size(audio.aligned_audioCube,2)
            %Warping
            x = 1:numel(brng);
            v = squeeze(dataBlock(j,k,:));
            tmp(j,k,:) = interp1(x,v,smp);
        end
        
    end
    
    %Load into carry-through vars
    ids = i.*ones(size(idx));
    sylType = [sylType, ids];
    sylSpecs = cat(1, sylSpecs, tmp);
    
    
    %Clear vars
    clear('audio', 'data')
end

%Check destination
if ~exist(dLoc, 'dir')
    mkdir(dLoc)
end

%Save the formatted data to file
save([dLoc, filesep, 'Bird XXX.mat'], 'sylType', 'sylSpecs', '-v7.3')








