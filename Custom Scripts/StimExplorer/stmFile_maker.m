%This script is designed to make a stimFile post-hoc (i.e., using a manual
%register of stim parameters and the folder of recordings during that
%time).

%Define the stimulation parameters that will be logged in a stimFile. We'll
%need to capture the start and end times/days and the stimulation configuration that
%corresponds to those periods. We're going to assume simple stimulation
%(i.e., common to all channels) for this stuff... 
formatIn = 'yyyy-mm-dd HH:MM:SS';
starts = {'2018-03-06 14:22:00';...
          '2018-02-26 08:00:00'; ...
          '2018-02-27 16:48:00'; ...
          '2018-03-05 10:25:00'; ...
          '2018-02-28 20:22:00'; ...
          '2018-03-01 11:35:00';...
          '2018-03-01 19:27:00'; ...
          '2018-03-02 17:18:00'};

ends =   {'2018-03-08 12:14:00';...
          '2018-02-27 16:48:00'; ...
          '2018-02-28 20:22:00'; ...
          '2018-03-06 14:21:00'; ...
          '2018-03-01 11:35:00'; ...
          '2018-03-01 19:27:00';...
          '2018-03-02 17:18:00'; ...
          '2018-03-05 10:25:00'};
      
I1 = [100, 100, 100, 100, 100, 100, 100, 100];
T1 = [100, 100, 100, 100, 100, 100, 100, 100];
TI = [5, 5, 5, 5, 5, 5, 5, 5];
I2 = [-100, -100, -100, -100, -100, -100, -100, -100];
T2 = [100, 100, 100, 100, 100, 100, 100, 100];
Reps = [5, 5, 5, 10, 10, 10, 20, 20];
Freq = [10, 50, 100, 10, 50, 100, 50, 100];

%Select the folders to use as "source data" -- will determine the
%"elements" entries to the stimFile
mother = 'V:\SongbirdData\LO58\';
folders = {'2018-02-26'; '2018-02-27'; '2018-02-28'; '2018-03-01'; '2018-03-02'; '2018-03-03'; ...
    '2018-03-04'; '2018-03-05'; '2018-03-06'};

numFolders = numel(folders);
for i = 1:numFolders
   %Get the folder listing
   listing = dir([mother, folders{i} '\*.stm']);

   %Simplify
   filename = getFieldVectorCell(listing', 'name');
   
    %Replicate the path to match the file listing length
paths = repmat([mother, folders{i}], numel(filename), 1);



end

checkIt = 1;


