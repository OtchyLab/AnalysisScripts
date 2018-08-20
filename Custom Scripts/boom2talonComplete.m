function boom2talonComplete
%function boom2talonComplete
%
%This function/script will reformat .WAV files ffrom GLab Boom Recorder to the
%Talon Standard (for compatibility with the rest of the processing
%pipeline). This version takes the raw files produced by the recording
%computer (~5.5MB), chops them using the zftftb scripts, renames files, and
%sorts them appropriately (i.e., by name and date)
%
% Created by TMO (1/6/2017);

%Set the locations for all file processing
rawLoc = 'C:\Users\Tim\Desktop\Song Data\Raw';
trashLoc = 'C:\Users\Tim\Desktop\Song Data\Trash';
finalLoc = 'C:\Users\Tim\Desktop\Song Data\Converted';

%Set the file conversion constants
fs = 44150; %desired final sample rate (Hz)
depth = 16; %bits per sample
refTime = datenum([1903 12 31 20 00 00]); %reference for timestamping
rmSourceFlag = false;

%Start processing loop

%Get the list of folders in the raw folder to be converted. Note that the
%foldernames found will become the bird names in later stages.
folderList = dir(rawLoc);
dirBool = getFieldVector(folderList', 'isdir');
dirCell = getFieldVectorCell(folderList(dirBool)', 'name');
folderList = dirCell(3:end); %parsed list of folder/birds to process

%Cycle trhough the folderList and sequentially convert to new format
for i = 1:numel(folderList)
    %Chop the long recordings in the current folder using the zftftb
    %scripts with the "trash" modificiation (i.e., output in trash
    %location)
    chopDir = [rawLoc filesep folderList{i}];
    zftftb_song_chop(chopDir)
    
    %Announce Step completion
    disp(['Completed chopping folder ' chopDir '. Now converting files...'])
    
    %If that process produced chopped files in the trash location, continue
    %processing and reformatting
    wavList = dir([trashLoc filesep '*.wav']);
    if ~isempty(wavList)
        for j = 1:numel(wavList)
            %Set current working file
            curFile = wavList(j).name;
            
            %Extract recording time from the filename
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
            
            WAVname = [folderList{i}, '_', serialNum, '_', fileYear, '_', fileMonth, '_', ...
                fileDay, '_', fileHour, '_', fileMin, '_', fileSec, '.wav'];

            %Identify desired file destination
            fileDate = [fileYear '-' fileMonth '-' fileDay];
            targetParent = [finalLoc filesep folderList{i}];
            targetFolder = [finalLoc filesep folderList{i} filesep fileDate];
            
            %If targetFolder doesn't exist, create it
            if ~exist(targetParent, 'dir')
                mkdir(targetParent)
            end
            
            if ~exist(targetFolder, 'dir')
                mkdir(targetFolder)
            end
            
            %Load wav from file
            readLoc = [trashLoc filesep curFile];
            [y, ~] = audioread(readLoc);

            %Write wav to destination with new filename
            writeLoc = [targetFolder filesep WAVname];
            audiowrite(writeLoc, y, fs, 'BitsPerSample', depth);
            
            %Delete the corresponding trash file
            delete(readLoc)
            
            %Announce Step completion
            disp(['Completed processing file ' WAVname])
        
        end
        
    end
    
    %If the flag is set, delete the source files from the raw location
    if rmSourceFlag
        rmdir(chopDir)
    end
    
end



















