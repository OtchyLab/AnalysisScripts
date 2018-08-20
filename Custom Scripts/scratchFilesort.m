%Function runs SAP similarity analysis on all wav-files in a folder. The batch file handling and actual similarity scoring is
%done by external scripts. Here I just handles loading from disk, formating, and calling subfunctions

%Select folders containing song snips to process
folders = uipickfiles('FilterSpec', 'V:\Rd279\', 'Prompt', 'Select folders to process');
numFolders = numel(folders);

mother = 'V:';

%Cycle through folders and extract out songs from each
for i = 1:numFolders
    %retrieve names of all wav-files in the folder
    files = dir([folders{i}, filesep, '*.wav']);
    
    %Retrieve foldername
    sp = regexp(folders{i}, '\', 'split');
    foldName = sp{end};
    
    %Sequentially load each file and parse get contens
    for j = 1:numel(files)
        %Get the bird name from file name
        sp = regexp(files(j).name, '_', 'split');
        birdName = sp{1};
        
        %Construct file paths
        source = [folders{i} filesep files(j).name];
        dest = [mother, filesep, birdName, filesep, foldName, filesep, files(j).name];
        
        %Only move if the location is different
        if ~strcmp(source, dest)
            %Create the derstination directory if it doesn't exist yet
            if ~exist([mother, filesep, birdName, filesep, foldName], 'dir')
                mkdir([mother, filesep, birdName, filesep, foldName])
            end
            
            movefile(source,dest)
        end
    end
    
    
end