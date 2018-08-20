function batchDAT6toDAT5
%This function is to batch reconstruct the original DAT-files from the
%DAT-files. This function handles the file indexing and user input, but all
%of the actual reconstruction is done by DAT6toDAT5.m

%Get the location of the original DAT-files to cross-reference
DATpath = uigetdir('V:\Pur238\','Select location of the source DAT files.');

%Generate the file listing for the DATs
DATfiles = dir([DATpath filesep '*.dat']);
if length(DATfiles) < 1
    disp(['There are no files to process in this folder. Check the location and try again.'])
    beep
    return
end

%Get the location to save the output DATs
%newDATpath = uigetdir('V:\Pur238\','Select location to save the reconstructed DAT files.');
newDATpath = DATpath; %This will overwrite the original files

%Loop through each DAT file, converting and saving as you go...
parfor i = 1:length(DATfiles)
    %Generate the full path of the DAT file
    fullFilename{i} = [DATpath filesep DATfiles(i).name];
    
    %Call function to recreate file and filename
    [filename, DAT{i}] = DAT6toDAT5(fullFilename{i});
    
    %Generate the full path of the new DAT file
    writeName{i} = [newDATpath filesep DATfiles(i).name];
    
    %Binary write
    fid{i} = fopen(writeName{i}, 'w');
    fwrite(fid{i}, DAT{i}, 'single', 0, 'b');
    fclose(fid{i});
    
    %Display write name to monitor progress
    disp([DATfiles(i).name ' written to file.'])
    
end




