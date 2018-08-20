%%% Lists all of the functions within a given folder; extracts ones who's
%%% names start with the name of the file
% Takes as an argument a cell array in which the first row contains
% foldernames and the second row contains the path names
function relevantFiles = getRelScripts(namesAndPaths)
    % change current directory to the pathname given
    relevantFiles = cell(1,1);
    count = 1;
    for j = 1:length(namesAndPaths)
        pathName = namesAndPaths{2,j};
        folderName = namesAndPaths{1,j}; 
        cd(pathName);
        desiredForm = strcat(folderName, "*.m");
        list = dir(desiredForm);
        %fprintf("%d\n%s\n", length(list), list(1).name)
        if length(list) > 1
            for i = 1:length(list)
                namelength = length(list(i).name);
                if namelength <= (strlength(folderName) +3)
                    relevantFiles{1,count} = list(i).name;
                    count = count + 1;
                end
            end
        else
            relevantFiles{1,count} = list.name;
            count = count + 1;
        end
        %fprintf("%s\n%s\n",folder_contents{1},folder_contents{2})

        % list files with names containing the name of the folder
    end
end