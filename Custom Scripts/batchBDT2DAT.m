function batchBDT2DAT
%This function is to batch reconstruct the original DAT-files from the
%BDT-files. This function handles the file indexing and user input, but all
%of the actual reconstruction is done by BDT2DAT.m

%Get the location of the original BDT-files to cross-reference
BDTpath = uigetdir('V:\Single Units\','Select location of the source BDT files.');

%Generate the file listing for the BDTs (restrict to only channel 0 to eliminate
%repeats)
BDTfiles = dir([BDTpath filesep '*0.BDT']);
if length(BDTfiles) < 1
    disp(['There are no files to process in this folder. Check the location and try again.'])
    beep
    return
end

%Get the location to save the output DATs
DATpath = uigetdir('V:\Single Units\','Select location to save the reconstructed DAT files.');

%Loop through each BDT file, converting and saving as you go...
for i = 1:length(BDTfiles)
    %Generate the full path of the BDT file
    fullFilename = [BDTpath filesep BDTfiles(i).name];
    
    %Call function to recreate file and filename
    [DATname,DAT] = BDT2DAT(fullFilename);
    
    %Generate the full path of the new DAT file
    writeName = [DATpath filesep DATname];
    
    %Binary write
    fid = fopen(writeName, 'w');
    fwrite(fid, DAT, 'single', 0, 'b');
    fclose(fid);
    
    %Display write name to monitor progress
    disp([DATname ' written to file.'])
    
end




