% folderSpecs.m
% 3/28/16 edited by TMO
%
% This function generates spectrograms for all .wav files found within specified folder. Standardized plotting functions are
% included. It assumes the standard .wav formating, but does not require the naming convention.
%

%Set the folder to scan
% folder = 'C:\Users\Tim\Desktop\Iso Paper Images\Example Tutored\045dph';
folder = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis\Example Spectrogram Sources';
% folder = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis';
% folder = '/Users/Tim/Documents/MATLAB/General/Custom Scripts/tCAF Analysis';

%Constants for Bandpass Audio (300-8000Hz)
fs = 44150;
HP_fNorm = 300/(fs/2);
LP_fNorm = 10000/(fs/2); %Changed 7/26/12 to match Farhan and Cengiz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
load('tCAF_cmap2.mat');

%Scan folder
files = dir([folder filesep '*.wav']);
% files = dir([folder filesep '*.stm']);
if isempty(files)
    print('Folder is empty. Check it or specify a new folder')
    return
end

%Cycle through the list to plot spectrograms
for i = 1:numel(files)
    %Load file from disk
    name = [folder filesep files(i).name];
    [raw, ~] = audioread(name);
%     d = getChannels(name);
%     raw = d(1,:);
    
    %Process the audio
    faudio = filtfilt(BP_b, BP_a, (raw-mean(raw))');
    
    %Plot spectrogram to a new figure
    figure('Name', files(i).name); clf
    displaySpecgramQuick(faudio, fs, [0,8000], [-9.5,12], 0);
    h = gca;
    h.Colormap = cmap;
    set(h, 'Box', 'off', 'TickDir', 'out')
    
    %Format the figure
    set(gcf, 'Units', 'Inches', 'Position', [0.5, 7.75, 13.5, 2.5])
    
end