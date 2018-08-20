%%%
% RUN THIS RIGHT AFTER FS_AV_Parse 


%%%%%%%%%%%%%%%%%%
% CSV 
%%%%%%%%%%%%%%%%%%

% Set working directory 
cd '/Users/sadiela/Desktop/goodVids'


% It works!!! :))))

% now do it in bulk... ugh 
% 
% % get list of all csv files
% csvfiles = dir(fullfile(bigDir, '*.csv'));
% [csvnames, index] = sortrows({csvfiles.name});
% 
% for i = 1:length(csvnames)
%     % get csv data in a cell array
%     data = csv2cell(csvnames{i}, 'fromfile');
%     % need to delete 7 rows starting from row 14 (14-20)
% 
%     data(14:24,:) = [];
%     % easy peasy! :))
%     [rows,cols] = size(data);
%     emptyCells = cellfun(@isempty,data);
% 
%     for j = 1:rows
%         for k = 1:cols
%             if emptyCells(j,k) == 1
%                 data{j,k} = "";
%             end
%         end
%     end
% 
% 
%     fid = fopen('2018-06-06 08 53 50.csv','w');
%     for l = 1:rows
%         fprintf(fid, '%s, %s, %s, %s\n', data{l,:});
%     end
%     fclose(fid); 
%     
%     clear rows cols fid data i j k l emptyCells
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%
% MAT
%%%%%%%%%%%%%%%%%%%%%%%%

% get list of mat files
matDir = '/Users/sadiela/Desktop/frameDeletion/mat';
bigDir = '/Users/sadiela/Desktop/frameDeletion';

matfiles = dir(fullfile(matDir, '*.mat'));
[matnames, index] = sortrows({matfiles.name}');

% load, edit, save, clear?? 
for i = 2:length(matnames)
    filePath = strcat(matDir, '/', matnames{i});
    load(filePath);
    % make video edits
    video.frames = video.frames(:,:,:,12:end);
    video.times = video.times(12:end);
    video.gain = video.gain(12:end);
    video.nrFramesTotal = length(video.gain);
    
    % save to same filename and clear variables
    save(filePath, 'audio', 'video'); 
    clear audio video filePath 
end

% STILL NEED TO CHANGE FILE NAMES!!!! (eventually)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = imread('2018-06-06 14 12 31_max_.png');






