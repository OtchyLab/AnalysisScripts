%rename files in the selected folder from *.dat to *.wav (nno other change to the file is made)

%Get the directories to analyze from user
folders = uipickfiles('FilterSpec', 'V:\Pur128\');

for i = 1:numel(folders)
    %Get all dats in the current folder
    files = dir([folders{i} filesep '*.wav']);
    
    %Cycle through files and rename
    for j = 1:numel(files)
        %old path
        oldName = [folders{i} filesep files(j).name];
    
        %new path just appends the wav suffix
        newName = [oldName(1:end-3) 'dat'];
        
        %move it
        movefile(oldName, newName);
        
    end
end