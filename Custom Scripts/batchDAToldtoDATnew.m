function batchDATold2DATnew
%This function takes as input the full filename (i.e., absolute path + filename) for an old-style DAT file (i.e.,
%sampled at 30,030Hz and not with the new name formatting) and reconstructs all four channels into the
%modern DAT format. The function returns both the reconstructed DAT and the updated filename. This function
%handles the file indexing and user input, but all of the actual reconstruction is done by DAT30k2DAT44k.m
%
%This will update the very oldest recordings TMO made, and those that BPO recorded in FeeLab

%Get the location of the original DAT-files to cross-reference
DATpath = uigetdir('V:\Pur127\','Select location of the source DAT files.');

%Generate the file listing for the DATs
DATfiles = dir([DATpath filesep '*.DAT']);
if length(DATfiles) < 1
    disp(['There are no files to process in this folder. Check the location and try again.'])
    beep
    return
end

%Get the location to save the output DATs
newDATpath = uigetdir('V:\Rd284\Neural Data\Fresh\','Select location to save the reconstructed DAT files.');
% newDATpath = DATpath; %This will overwrite the original files

%Loop through each DAT file, converting and saving as you go...
for i = 1:length(DATfiles)
    %Generate the full path of the DAT file
    fullFilename = [DATpath filesep DATfiles(i).name];
    
    %Call function to recreate file and filename
    [DATname, DAT] = DATold2DATnew(fullFilename);
    
    %Generate the full path of the new DAT file
    writeName = [newDATpath filesep DATname];
    
    %Binary write
    fid = fopen(writeName, 'w');
    fwrite(fid, DAT, 'single', 0, 'b');
    fclose(fid);
    
    %Display write name to monitor progress
    disp([DATfiles(i).name ' written to file.'])
    
end




