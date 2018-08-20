function [files]=FilesList
%Simply returns a variable that contains the names of the Feature files in
%the current folder

Listed = dir('*features.mat');

for i=1:length(Listed)
    buffer{i} = Listed(i).name;
end

files = buffer;
end