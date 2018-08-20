function [filename, DAT] = DAT6toDAT5(fullFilename)
%This function takes as input the full filename (i.e., absolute path +
%filename) for a 6-channel DAT file and convertas it to the standard 5-channel DAT structure.
%The function returns the reconstructed DAT.

%For simplicity, this is the test filename; comment out later.
%fullFilename = 'V:\Pur238\2010-05-25\Pur238_3357652217_2010_05_25_13_10_17.dat';

% %Parse full filename using the fileseps, and extract filename
splitFolders = regexp(fullFilename, filesep, 'split');
filename = splitFolders{end};

%Preallocate variables
chunks = [];
channels = [];

%Read the file
fid = fopen(fullFilename);

%Parse header information
fs = fread(fid, 1,'single',0,'b');
fs = 44150; %Hardcode this for now as many of the headers are incorrect.
crap = fread(fid, 1,'single',0,'b');
numChan = fread(fid, 1,'single',0,'b');
numChan = 5;

%The first one second of data from each channel is read sequentially
for i=1:numChan
   chunks(i,:) =  fread(fid, fs,'single',0,'b');
end

%The the data is interleaved to the end of the file
MultiPlexed = fread(fid, inf,'single',0,'b');
for i=1:numChan
   %data(i,:) = [chunks(i,:), MultiPlexed(i:numChan:end)'];
   data(i,:) = [chunks(i,:), MultiPlexed(i:6:end)'];  %Correction for 6-chan mistake
end
fid = fclose('all');

%Use the original decoding scheme to reassemble into a single dim vector
chunkSize = 44150; %The correct chunk size is 44150, but the original BDT encoding was fucked up in some cases. fix it with 44100
[m,n] = size(data);
endSeg = zeros(1,(m*n)-(m*chunkSize));
initSeg = [data(1,1:chunkSize), data(2,1:chunkSize), data(3,1:chunkSize), data(4,1:chunkSize), data(5,1:chunkSize)];
for i = 1:5
    endSeg(i:5:end) = data(i,(chunkSize+1):end);
end

%Add the header information to the final array
DAT = [44150, 420, 5, initSeg, endSeg];

%To write this to file properly, use:
%     fid = fopen(filename, 'w');
%     fwrite(fid, DAT, 'single', 0, 'b');
%     fclose(fid);
