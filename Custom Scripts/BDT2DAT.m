function [DATname,DAT] = BDT2DAT(fullFilename)
%This function takes as input the fill filename (i.e., absolute path +
%filename) for an old-style BDT file and reconstructs all four channels into the
%original DAT format. The function returns both the reconstructed DAT and
%the updated filename.

%For simplicity, this is the test filename; comment out later.
%fullFilename = 'V:\Single Units\Pur172 Sorted Cells\Pur172 091010 392dph\Pur172_d000038_20100910T131004chan0.BDT';

%Parse full filename using the fileseps, and extract filename
splitFolders = regexp(fullFilename, filesep, 'split');
filename = splitFolders{end};

%Parse filename using the common naming structure
splitName = regexp(filename,'_','split');
birdName = splitName{1};
dateStr = splitName{3};
fileYear = (dateStr(1:4));
fileMonth = (dateStr(5:6));
fileDay = (dateStr(7:8));

fileHour = (dateStr(10:11));
fileMin = (dateStr(12:13));
fileSec = (dateStr(14:15));

%Calculate the unique timestamp/serial number
refTime = datenum([1903 12 31 20 00 00]);
timeStamp = datenum(str2num(fileYear), str2num(fileMonth), str2num(fileDay),...
    str2num(fileHour), str2num(fileMin), str2num(fileSec));
serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);

%Reconstruct the name from components
DATname = [birdName, '_', serialNum, '_', fileYear, '_', fileMonth, '_', ...
    fileDay, '_', fileHour, '_', fileMin, '_', fileSec, '.dat'];

%Based on the passed filename, load from file all five channels of data
%into a single matrix
baseName = fullFilename(1:end-5);
for i = 1:5
    [~, data(i,:), ~, ~]  = daq_readDatafileBence([baseName num2str(i-1) '.BDT']);
end

%Use the original decoding scheme to reassemble into a single dim vector
chunkSize = 44100; %The correct chunk size is 44150, but the original BDT encoding was fucked up in some cases. fix it with 44100
[m,n] = size(data);
endSeg = zeros(1,(m*n)-(m*chunkSize));
initSeg = [data(1,1:chunkSize), data(2,1:chunkSize), data(3,1:chunkSize), data(4,1:chunkSize), data(5,1:chunkSize)];
% initSeg0 = [data(1,1:chunkSize), data(2,1:offset)];
% initSeg1 = [data(2,offset:44100) data(3,(2*offset):44100)];
% initSeg2 = [data(3,(3*offset):44100) data(4,(2*offset):44100)];
% initSeg3 = [data(2,offset:44100) data(3,(2*offset):44100)];
for i = 1:5
    endSeg(i:5:end) = data(i,(chunkSize+1):end);
end

%Add the header information to the final array
DAT = [44150, 420, 5, initSeg, endSeg];

%To write this to file properly, use:
    %fid = fopen(DATname, 'w');
    %fwrite(fid, DAT, 'single', 0, 'b');
    %fclose(fid);
