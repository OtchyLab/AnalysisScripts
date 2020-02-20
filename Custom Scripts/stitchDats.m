%DAT file stitching script
%
% Created by TMO

clear all

%Files to stitch
mother = '/Users/Tim/Dropbox/Test Recordings';
f{1} = 'LR81RY177_3586769480_2017_08_28_08_51_20.402.dat';
%f{2} = 'LR81RY177_3586464005_2017_08_24_20_00_05.008.dat';
f{2} = 'LR81RY177_3586465057_2017_08_24_20_17_36.511.dat';
f_out = 'LR81RY177_3586769489_2017_08_28_08_51_29.612.dat';

%Sequential load
data = [];
for i = 1:numel(f)
    %Full path name
    file = [mother, filesep, f{i}];
    
    %Load from file
    [chans,fs] = getChannels(file);
    
    %Stitch
    if i==1
        data = [data, chans(:,1.062e4:3.184e5)];
    else
        data = [data, chans];
    end
end

%Use the original decoding scheme to reassemble into a single dim vector
chunkSize = fs;
[m,n] = size(data);
endSeg = zeros(1,(m*n)-(m*chunkSize));
initSeg = [data(1,1:chunkSize), data(2,1:chunkSize), data(3,1:chunkSize), data(4,1:chunkSize), data(5,1:chunkSize)];
for i = 1:5
    endSeg(i:5:end) = data(i,(chunkSize+1):end);
end

%Add the header information to the final array
DAT = [44150, 420, 5, initSeg, endSeg];

%Write back to file
writeName = [mother, filesep,  f_out];

%Binary write
fid = fopen(writeName, 'w');
fwrite(fid, DAT, 'single', 0, 'b');
fclose(fid);