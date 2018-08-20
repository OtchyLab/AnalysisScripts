% get the path of a folder in a specified directory
function namePathCell = getPathNames(foldNames, baseDir, opt)
% option 1: creates a 1xlength(foldNames) cell array containing the path
% names for all folders
% option 2: creates a cell array with two rows and length(foldern) columns.
% The first row contains the folder names and the second row contains the path
% names. 
if opt == 1
    namePathCell = cell(1,length(foldNames));
        for i = 1:length(foldNames)
            %name_path_cell{1,i} = foldnames{i};
            namePathCell{1,i} = strcat(baseDir, "/", foldNames{i});
        end
elseif opt == 2
    namePathCell = cell(2,length(foldNames));
        for i = 1:length(foldNames)
            namePathCell{1,i} = foldNames{i};
            namePathCell{2,i} = strcat(baseDir, "/", foldNames{i});
        end
end
    
end