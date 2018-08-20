function boom2talon
%function boom2talon
%
%This function/script will reformat .WAV files ffrom GLab Boom Recorder to the
%Talon Standard (for compatibility with the rest of the processing
%pipeline). 
%
% Created by TMO (11/2/2016); last modified by TMO (11/3/2016)

%Set the file conversion constants
fs = 44150; %sample rate (Hz)
depth = 16; %bits per sample
refTime = datenum([1903 12 31 20 00 00]); %reference for timestamping

%Specify the location of files to convert
birdname = 'llb59rblk132';
source = 'C:\Users\Tim\Desktop\Song Data\';
dest = 'C:\Users\Tim\Desktop\Song Data\Converted\';

%Create destination folder if it does not exist
if ~isdir([dest birdname])
    mkdir([dest birdname])
end

%Capture the directory structure to be converted
dateFoldersStruct = dir([source birdname filesep]);
dirBool = getFieldVector(dateFoldersStruct', 'isdir');
dirCell = getFieldVectorCell(dateFoldersStruct(dirBool)', 'name');
dirCell = dirCell(3:end); %remove garbage

%Cycle through the directories and convert sequetially
for i = 1:numel(dirCell)
% for i = 1
    %Confirm it's of the right structure and not empty
    structChk = isdir([source birdname filesep char(dirCell(i)) filesep 'wav' filesep]);
    if structChk
        wavList = dir([source birdname filesep char(dirCell(i)) filesep 'wav' filesep 'c*.wav']);
        stuffChk = ~isempty(wavList);
    else
        stuffChk = false;
    end
    
    %If you meet the prereqs, convert the files in the folder
    if structChk && stuffChk
        %Create destination folder if it does not exist
        if ~isdir([dest birdname filesep char(dirCell(i))])
            mkdir([dest birdname filesep char(dirCell(i))])
        end

        for j = 1:numel(wavList)
            %Set current file
            curFile = wavList(j).name;
            
            %Extract the recording time from the filename
            splitName = regexp(curFile,'_','split');
            dateStr = splitName{2};
            fileYear = dateStr(4:7);
            fileMonth = dateStr(9:10);
            fileDay = dateStr(12:13);
            fileHour = dateStr(15:16);
            fileMin = dateStr(18:19);
            
            bite = str2num(dateStr(end)); chunk = str2num(splitName{4}(1));
            offset = bite*40 + chunk*6;
            fileMin = num2str(str2num(fileMin) + floor(offset/60), '%02d');
            fileSec = num2str(round(60*((offset/60) - floor(offset/60))), '%02d');
            
            %Generate a Talon-compatible filename from the recording time
            timeStamp = datenum(str2num(fileYear), str2num(fileMonth), str2num(fileDay),...
                str2num(fileHour), str2num(fileMin), str2num(fileSec));
            serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);
            
            WAVname = [birdname, '_', serialNum, '_', fileYear, '_', fileMonth, '_', ...
                fileDay, '_', fileHour, '_', fileMin, '_', fileSec, '.wav'];

            %Load wav from file
            readLoc = [source birdname filesep char(dirCell(i)) filesep 'wav' filesep curFile];
            [y,Fs] = audioread(readLoc);
            
            %Write wav to destination with new filename
            writeLoc = [dest birdname filesep char(dirCell(i)) filesep WAVname];
            audiowrite(writeLoc, y, fs, 'BitsPerSample', depth);
        end
        
        %Announce Step completion
        display(['Completed processing folder ' source birdname filesep char(dirCell(i))])
    else
        display(['Skipping folder ' source birdname filesep char(dirCell(i))])
    end
end


