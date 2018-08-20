% function to extract dependencies of a particular file and store them in a
% new folder
%%% arguments: cell array of file names, cell array of folder paths, target
%%% directory
% this function contains several subfunctions only used by it
function generateDependenceList(fileNames, folderPaths, targetPath)
    % set working directory to '/Users/labperson/Desktop/Software'
    folderPaths{1,length(folderPaths)+1} = targetPath;
    test_paths = folderPaths;
    cd('/Users/labperson/Desktop/Software');
    for k = 1:length(fileNames)
    dep = getDependencies(fileNames{k});
    %displayDependencies(dep)
    % copy into target directory
        for i = 1:length(dep)
            % CHECK IF FILE IS ALREADY IN TARGET DIRECTORY!!!
            % Check if file is in a desired package directory
            if ~checkfolders(dep{1,i}, test_paths)
                copyfile(dep{1,i}, targetPath) % this works!!! :))
                fprintf("%s\n",dep{1,i});
            end
        end
    end
end
%%%%%%%%%%%% subfunctions %%%%%%%%%%%%
% prints filepaths of all dependencies 
function displayDependencies(dependencies)
    for i = 1:length(dependencies)
        fprintf("%s\n",dependencies{1,i});
    end
end
% gets paths of all dependencies of a scripts and puts them in a cell array
function dependencies = getDependencies(fileName)
% returns pathnames of all dependencies of a given file
    [dependencies, pList] = matlab.codetools.requiredFilesAndProducts(fileName);
end
function contains = checkfolders(filePathName, testPathNames)
% returns true if the common files folder contains a file with the specified
% name at the end of the path input OR if the file is in one of the source
% folders
% arguments: path of file in question, cell array of directories of which 
% to test membership
    contains = 0;
    count = 1;
    while contains == 0 && count <= length(testPathNames)
        fileName = getFileName(filePathName);
        %empty if file does not yet exist (false)
        %fprintf("%s\n%s\n", test_pathnames{count}, filename)
        dirstruct = dir(fullfile(testPathNames{count}, fileName));
        if length(dirstruct) == 1
            contains = 1;
        end
        count = count +1;
    end
end
function fname = getFileName(pathName)
% returns a file name when given a full path name
    i = length(pathName);
    while ~strcmp(pathName(i), "/")
        i = i-1;
    end
    fname = pathName(i+1:length(pathName));
end