function name = getBirdName(path)
    %filepath = fileparts(which(fileName));
    endIndex = length(path);
    startIndex = endIndex; 
    while (~strcmp(path(startIndex), '/'))
        startIndex = startIndex - 1; 
    end
    name = path(startIndex+1:endIndex);
end