function [cellVect,doubleVect] = getKeysSubset(keys,field)
% Assumes that the input 'keys' is an array of strings specifying the
% files names for all the records in the annotation file.  'Keys' has
% the form: 'BirdName_SerialNum_YYYY_MM_DD_HH_MM_SS.FileType'.  Valid codes
% for 'field' are:
% 'name' -- BirdName
% 'serial' -- SerialNum
% 'y' -- year
% 'mon' -- month
% 'd' -- day
% 'h' -- hour
% 'min' -- minute
% 's' -- second
% 'type' -- FileType

% Decode 'field' to a pointer value
if strcmp(field,'name')
    pntr = 1;
elseif strcmp(field,'serial')
    pntr = 2;
elseif strcmp(field,'y')
    pntr = 3;
elseif strcmp(field,'mon')
    pntr = 4;
elseif strcmp(field,'d')
    pntr = 5;
elseif strcmp(field,'h')
    pntr = 6;
elseif strcmp(field,'min')
    pntr = 7;
elseif strcmp(field,'s')
    pntr = 8;
elseif strcmp(field,'type')
    pntr = 9;
else
    error('Invalid field specified -- check syntax.')
    return
end

%Initialize output variables
cellVect = [];
doubleVect = [];

%Loop through the whole list to extract the specified field
numStr = length(keys);
for i = 1:numStr
    %Parse file name for components
    temp = regexp(keys{i},'_','split');
    temp{9} = temp{8}(end-2:end);
    temp{8} = temp{8}(1:end-4);
    
    %Populate output vectors
    cellVect{i} = temp{pntr};
    if pntr~=1 || pntr~=9
       doubleVect(i) = str2double(temp{pntr});
    end
end

