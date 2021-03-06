% This script removes from an annotation file all of the old records that point to deleted/non-existant files.
%
% All input/output parameters are hard coded below.
% 
% TMO, 12/16/2014
% 
% Load Data

% Set the location of the input and output files/folders
annotName = 'X:\Sil167\Sil167_100407_annotation2.mat';
outputName = 'X:\Sil167\Sil167_100407_annotation2.mat';

dataDir = 'X:\Sil167\2010-04-07';

% Clean load annotation
elements = []; keys = [];
load(annotName);

% Generate directory file list
list = dir([dataDir, '*.*']);
files = cell(length(list),1);
for i = 1:length(list)
    files{i} = list(i).name;
end

% Generate indices of present files
indx = []; pntr = [];
for i = 1:length(keys)
    if ismember(keys{i}, files)
        indx(i) = true;
        pntr(end+1) = i;
    else
        indx(i) = false;
    end
end

% Only retain the records for files that exist
elements = elements(1,pntr);
keys = keys(1,pntr);

% Save at the new location
save(outputName, 'elements', 'keys', '-v7')

clear all;

disp('All done!')


