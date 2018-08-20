function [DATname,DAT] = DATold2DATnew(fullFilename)
%This function takes as input the full filename (i.e., absolute path + filename) for an old-named DAT file
%and reconstructs all four channels into the modern DAT format. The function returns both the
%reconstructed DAT and the updated filename.
%
%This will update the very oldest recordings TMO made, and those that BPO recorded in FeeLab

%For simplicity, this is the test filename; comment out later.
% fullFilename = 'V:\Rd284\Neural Data\62\Rd284_090730_133053.DAT';

%Parse full filename using the fileseps, and extract filename
splitFolders = regexp(fullFilename, filesep, 'split');
filename = splitFolders{end};

%Parse filename using the common naming structure
splitName = regexp(filename,'_','split');
birdName = splitName{1};
dateStr = splitName{2};
fileYear = ['20', (dateStr(1:2))];
fileMonth = (dateStr(3:4));
fileDay = (dateStr(5:6));

timeStr = splitName{3};
fileHour = (timeStr(1:2));
fileMin = (timeStr(3:4));
fileSec = (timeStr(5:6));

%Calculate the unique timestamp/serial number
refTime = datenum([1903 12 31 20 00 00]);
timeStamp = datenum(str2num(fileYear), str2num(fileMonth), str2num(fileDay),...
    str2num(fileHour), str2num(fileMin), str2num(fileSec));
serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);

%Reconstruct the name from components
DATname = [birdName, '_', serialNum, '_', fileYear, '_', fileMonth, '_', ...
    fileDay, '_', fileHour, '_', fileMin, '_', fileSec, '.dat'];

%Based on the passed filename, load from file all five channels of data into a single matrix
%Preallocate variables
chunks = [];
channels = [];

%Read the file
fid = fopen(fullFilename);
%dataset = fread(fid, inf,'single',0,'b');
%fclose(fid);

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
   channels(i,:) = [chunks(i,:), MultiPlexed(i:numChan:end)'];
   %channels(i,:) = [chunks(i,:), MultiPlexed(i:6:end)'];  %Correction for 6-chan mistake
end
fid = fclose('all');

%
data = channels;

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
    %fid = fopen(DATname, 'w');
    %fwrite(fid, DAT, 'single', 0, 'b');
    %fclose(fid);
