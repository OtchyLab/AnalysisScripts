function AudioStripper(type,loc)
%This function strips strips the audio out of DAT files and saves it as a
%WAV file in the same folder with the same filename + WAV suffix.  
%
%The input options are:
%   type: a string of either 'multidirs', 'folder' or 'file'.
%   location: if present, will set where the ui's start from

%Preset the input arguments
if nargin == 1
    loc = [];
end

if nargin == 0
    loc = [];
    type = 'file';
end

%Preallocate the vars
fs = 44150;
depth = 16;
folders.isdir = false;

%Get the filenames for the DATs to be processed
cd('V:\')
if strcmp(type,'file')
    %Process one or more DAT files by name
    [temp,path] = uigetfile([loc '*.dat'],'Select the files to convert','MultiSelect','on');
    cd(path);
    if size(temp,2)>1
        for i = 1:size(temp,2)
            fnames(i).name = char(temp{i});
        end
    else
        fnames.name = temp;
    end
    children = 1;
elseif strcmp(type,'folder')
    %Process all DAT files in a folder
    path = uigetdir(loc,'Select the folder to convert');
    cd(path);
    fnames = dir('*.dat');
    children = 1;
elseif strcmp(type,'multidirs')
    %Process all DAT files across several directories
    parent = uigetdir(loc,'Select the parent folder to start at');
    cd(parent)
    folders = dir('20*'); %Assuming it's the usual date naming structure
    children = length(folders);
end

%If in multi-mode, step through each child folder
for i = 1:children
    if strcmp(type,'multidirs') && folders(i).isdir
        %get the filenames for all the DATs in the child folder
        cd([parent '\' folders(i).name])
        fnames = dir('*.dat');
    elseif strcmp(type,'multidirs') && ~folders(i).isdir
        fnames = [];
    end
    
    %Work through each file in the fname structure
    for j = 1:length(fnames)
        %Load the file into memory and strip the channels out of
        %the raw recording file
        [rawdata, ~] = getChannels(fnames(j).name);
        audio = rawdata(1,:);
        %figure(1)
        %plot(audio)

        %Check for interesting audio content: syllables vs noise/crap)
%        [syllStartTimes, ~, ~, ~, ~, ~, ~, ~] = aSAP_segSyllablesFromRawAudioBence(audio', fs, 0.5, 0.7);
%        if length(syllStartTimes) >= 3
            convert = true;
%        else
%            convert = false;
%        end

        %Write the audio data to the new wave file
        if convert
            newname = [fnames(j).name(1:end-3) 'wav'];
            wavwrite(audio,fs,depth,newname);
        end

        %Delete the original dat file
        %delete([fnames(j).name]);
    end
    fnames = [];
end
