function remapFolders
%This function grabs all files (meeting a criterion) in the main directory and all subdirectories. It calls a function
%(at this point a subfunction to this file) that generates a new files name based on the old, and moves it to that location.

%Because the files from Ofer were so flooded with line noise as to be useless in the current state, I do some basic filtering
%(zero-phase Butterworth bandpass, 300-10000Hz) before saving them to the destination.

%3/17/2015 -- TMO

%Place to start search
topLevel = 'F:\day 43 samba\photoperiod 16 to 8\R658';
birdNum = '658';
%Set output location (later, just to V:\workspace)
outputLoc = 'C:\Users\Tim\Desktop\TempDrop';

%Generate conversion log file in the home destination directory
fileID = fopen([outputLoc filesep 'Bird' birdNum ' convertLog.txt'],'w');
fprintf(fileID,'%1s %50s %75s\r\n','Age','Old Name', 'New Name');

%Filtering constants (300-10000Hz)
HP_fNorm = 300/(44150/2);
LP_fNorm = 10000/(44150/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);

%Determine the structure of the underlying directory to traverse
topDir = dir(topLevel);
topDir = topDir(3:end); %remove '.' and '..'

dirIdx = logical(getStructField(topDir,'isdir'));
topDir = topDir(dirIdx); %Only leave the folders (I don't think there are stray files, but this is a good check)

%Loop through each of the directories and process sequentially
numFolders = length(topDir);
for i = 1:numFolders
    curDir = [topLevel, filesep, topDir(i).name];
    
    %Retrieve the wav-filelist from the current directory
    curList = dir([curDir filesep '*.wav']);
    

    %Generate new filename based on the old one and the standard nameing nomenclature
    targetLoc = [];
    fileYear = year(curList(1).date);
    for j = 1:length(curList)
        %Generate new filename
        newName = remapTch(curList(j).name, birdNum, fileYear);   
    
        %If this is the first file in the folder, determine the target location (directory)
        if j == 1 || isempty(targetLoc)
            %Use the filename to generate the folder name
            sp = regexp(newName, '_', 'split');
            folderDest = [char(sp{1}) filesep char(sp{3}), '-' char(sp{4}), '-' char(sp{5})];
            
            %Create the directory is necessary
            if exist([outputLoc filesep folderDest], 'dir') ~= 2
                mkdir([outputLoc filesep folderDest]);
            end
        end
        
        %Assemble the to/from locations for easy code reading
        curLoc = [curDir filesep curList(j).name];
        targetLoc = [outputLoc filesep folderDest filesep newName];
        
        %Load file from disk
        [audio, ~] = audioread(curLoc);
        
        %Bandpass filter
        audioFilt = filtfilt(BP_b,BP_a,audio);
        
        %Write filtered version to destination
        audiowrite(targetLoc, audioFilt, 44150);
        
        %Write to the text log file what exactly was done
        fprintf(fileID, '%1s %50s %70s\r\n', topDir(i).name, curList(j).name, newName);
        
    end
end

%Close the open file
fclose(fileID);

function newName = remapTch(oldName, birdNum, fileYear)
%Subfunction handles the parsing of the old file name and the generation of a new one.

%Prefix that will substitute for the band color inthe typical naming scheme.
prefix = 'Tch';

%Passed from calling program b/c these don't specify
fileYear = num2str(fileYear);
birdNum = birdNum;

%Split up and start parsing
sp = regexp(oldName(1:end-4), '_', 'split');

%Generate new bird name

birdName = [prefix, birdNum];

if strcmp(sp{3}, 'on') %Flag for the worded date format
    %Straight forward Copy
    fileDay = char(sp{5});
    fileHour = char(sp{6});
    fileMin = char(sp{7});
    
    %Copy or make it up
    if length(sp) > 7
        fileSec = char(sp{8});
    else
        fileSec = num2str(randi(59)); %Randomly assign seconds value to assure non-duplication for those files that don't specify a seconds field
    end
    
    %Use lookup table function
    fileMonth = num2str(name2num(sp{4}));

else
    %Straight forward Copy
    fileMonth = char(sp{3});
    fileDay = char(sp{4});
    fileHour = char(sp{5});
    fileMin = char(sp{6});
    
    %Copy or make it up
    if length(sp) > 6
        fileSec = char(sp{7});
    else
        fileSec = num2str(randi(59)); %Randomly assign seconds value to assure non-duplication for those files that don't specify a seconds field
    end
    
end

%Now that all the time/date info is gathered, calculate the unique timestamp/serial number
refTime = datenum([1903 12 31 20 00 00]);
timeStamp = datenum(str2num(fileYear), str2num(fileMonth), str2num(fileDay),...
    str2num(fileHour), str2num(fileMin), str2num(fileSec));
serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);

%Reconstruct the name from components
newName = [birdName, '_', serialNum, '_', fileYear, '_', fileMonth, '_', ...
    fileDay, '_', fileHour, '_', fileMin, '_', fileSec, '.wav'];


function monNum = name2num(monName)
%Lookup table turning name of a month to a number
names = [{'January'}, {'February'}, {'March'}, {'April'}, {'May'}, {'June'}, {'July'}, {'August'}, {'September'}, {'October'}, {'November'}, {'December'}];
nums = 1:12;

indx = strcmp(monName,names);

monNum = nums(indx);


























