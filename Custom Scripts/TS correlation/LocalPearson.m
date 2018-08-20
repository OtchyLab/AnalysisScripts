function rho = LocalPearson(TS1,TS2,winsizeT)

if length(TS1)~=length(TS2)
    print('The two time series must have the same length')
    return
end

% Windowing parameters
fs = 1000;%44150; % Sampling frequency in Hz
%winsizeT = 50; % Window size in ms
winstepT = 1; % Window step in ms

winsize = floor(winsizeT*(fs/1000)); % Convert to sample counts
winstep = floor(winstepT*(fs/1000)); % Convert to sample counts

% Window the two time series
Mat1 = running_windows(TS1,winsize,winstep);
Mat2 = running_windows(TS2,winsize,winstep);

winnum = size(Mat1,1);

for i = 1:winnum
    
    rho(i) = (corr(Mat1(i,:)',Mat2(i,:)','type','Pearson'));
end

for j = 1:floor((winsize/2))
    Bstart = max(1,j-26);
    Bend = j+24;
    segB = Bstart:Bend;
    rhoB(j) = (corr(TS1(segB)',TS2(segB)','type','Pearson'));
    
    Estart = length(TS1)-26+j;
    Eend = length(TS1);
    segE = Estart:Eend;
    rhoE(j) = (corr(TS1(segE)',TS2(segE)','type','Pearson'));
end

rho = [rhoB,rho,smooth(fliplr(rhoE(2:end)),5)'];

    