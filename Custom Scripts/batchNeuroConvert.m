function batchNeuroConvert(type, loc)
%This function batch processes neural recording files and converts them to the modern structure (i.e., 5 channels,
%44150 Hz sampling, unified naming convention). This top level function handles processing multiple folders of data
%and several types of conversion.
%
%Each conversion is handles by a separat subfunction, so changes tot he conversion should be handled at that level.
%
% ----Inputs----
%   type            This specifies the type of conversion to execute. This is a string input and valid inputs are:
%                     '6to5'     -- Converts fucked up 6-channel recordings to 5-channels; naming unchanged
%                     '30to44' -- Converts 30.03 recordings to 44.15kHz; also updates naming to unified
%                     'BDTtoDAT' -- Reassembles the 5 channels of BDT files to the original DAT structure
%                     'oldtonew' -- No changes to the data structure, simply renames to the unified format
%                     'doubleChunk' -- Fixes old error in which the first 1sec chunk of data was encoded as double instead of single precision
%                     '8to5'     -- Converts fucked up 8-channel recordings to 5-channels; naming unchanged
%
%   loc            This specifies where to save converted files. This is a string input and valid inputs are:
%                     'same'     -- Saves files to the same location. If there is no name change, old will be over-written
%                     'auto' -- Determines where to save the data based on the unified naming structure
%                     'ask' -- For each folder you process, you will be asked where you'd like to save Dats to
%
%This will update the very oldest recordings TMO made, and those that BPO recorded in FeeLab
%
%Last updated 6/24/2015 TMO

%Set flags for the input arguments
if nargin < 2
    disp('Missing required input arguments. Include them in the command and try again.')
    return
end

if strcmp(type, '6to5')
    conv = 1;
elseif strcmp(type, '30to44')
    conv = 2;
elseif strcmp(type, 'BDTtoDAT')
    conv = 3;
elseif strcmp(type, 'oldtonew')
    conv = 4;
elseif strcmp(type, 'doubleChunk')
    conv = 5;
elseif strcmp(type, '8to5')
    conv = 6;
else
    disp('No idea what type of conversion you want... check the "type" argument and try again.')
    return
end

if strcmp(loc,'same')
    fileLoc = 1;
elseif strcmp(loc, 'auto')
    fileLoc = 2;
elseif strcmp(loc,'ask')
    fileLoc = 3;
else
    disp('Where the hell am I saving this stuff? Check the "loc" argument and try again.')
    return
end

%Get the location of the original DAT-files to convert
DATfolders = uipickfiles('FilterSpec', 'V:\', 'Prompt', 'Select the folders with DATs to convert.');

if iscell(DATfolders)
    numFolders = length(DATfolders);
else
    disp('Something weird with the folder list... check that shit out, yo.')
end

for i = 1:numFolders
    %Generate the file listing for conversion
    if conv ~= 3
        DATfiles = dir([DATfolders{i} filesep '*.stm']); %change here...
    else 
        DATfiles = dir([DATfolders{i} filesep '*0.BDT']);
    end
    
    if length(DATfiles) < 1
        %If there are no files, report it and moveon
        disp(['There are no files to process in folder: ' DATfolders{i} '. Skipping now, but check the location.'])
        beep
    else
        
        %If "ask" or "same" location was selected, only do it at the start of processing the folder
        if fileLoc == 1
            newDATpath = DATfolders{i};
        elseif fileLoc == 3
            newDATpath = uigetdir('V:\',['Select location to save the reconstructed DAT files for folder: ' DATfolders{i} '.']);
        end
        
        %Loop through each DAT file, converting and saving as you go...
        h = waitbar(0,['Processing all files in folder: ' DATfolders{i} '. Please wait...']);
        for j = 1:length(DATfiles)
            %Generate the full path of the DAT file
            fullFilename = [DATfolders{i} filesep DATfiles(j).name];
            
            %Call relevant function to convert file
            if conv == 1
                [DATname, DAT] = DAT6toDAT5(fullFilename);
            elseif conv == 2
                [DATname, DAT] = DAT30k2DAT44k(fullFilename);
            elseif conv == 3
                [DATname, DAT] = BDT2DAT(fullFilename);
            elseif conv == 4
                [DATname, DAT] = DATold2DATnew(fullFilename);
            elseif conv == 5
                [DATname, DAT] = double2singleDAT(fullFilename);
            elseif conv == 6
                [DATname, DAT] = DAT8toDAT5(fullFilename);
            end

            %If "auto" location is selected, determine it based on the filename and the conventional file structure
            if fileLoc == 2 && j ==1
                sp = regexp(DATname,'_', 'split');
                newDATpath = ['V:\' char(sp{1}) filesep char(sp{3}) '-' char(sp{4}) '-' char(sp{5})];
                
                %Create folder if necessary
                if exist(newDATpath) ~= 7
                    mkdir(newDATpath);
                end
            end

            %Assemble the full filename to write
            writeName = [newDATpath filesep DATname];
            
            %Binary write
            fid = fopen(writeName, 'w');
            fwrite(fid, DAT, 'single', 0, 'b');
            fclose(fid);
            
            %Update waitbar
            waitbar(j/length(DATfiles))

        end
        
        %Close out the waitbar when the folder is done
        close(h)
    end
end










