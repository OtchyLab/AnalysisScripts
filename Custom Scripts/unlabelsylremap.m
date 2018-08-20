function [remapedelements] = unlabelsylremap(elements)
%This function remaps the randitional '-1' unknown segment type to #8 so
%that the files can be processed by the showRasters program more easily.

numFiles = length(elements);

for i=1:numFiles
    newType = elements{i}.segType;
    unkIndex = find(newType==-1);
    newType(unkIndex) = 8;
    elements{i}.segType = newType;
end

remapedelements = elements;

end