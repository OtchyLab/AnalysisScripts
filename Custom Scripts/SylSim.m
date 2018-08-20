function [SylParam, z_Params, raw_Params]=SylSim(songsource, dest)
%The most top-level function call for the syllable similarity program

%Ask User where to get source data
%songsource = uigetdir(pwd,'Source Data');
%songsource = 'C:\Documents and Settings\Bird Brain\Desktop\Test\Wave Chopping Test Folder\Source Data\';
cd(songsource);
files = dir('*.wav');

%Ask user where to store syllables
%dest = uigetdir(pwd,'Destination');
%dest='C:\Documents and Settings\Bird Brain\Desktop\Test\Wave Chopping Test Folder\Destination\';

%Batch call for chopping waves into syllables
for k = 1:numel(files);
    hour = regexp(files(k).name, '_', 'split'); %Parse file name for hour
    hour=str2double(hour(6));
    if hour>=12 && hour<=14; %Only process song between noon and 2:59pm
        file = files(k).name;
        path = [songsource '\' file];
        WaveChopper(path, dest);
    end
end

%Change working directory to the one containing the syllables and list file
%names
cd(dest);
files = dir('*.wav');

%Initialize variables; preallocate where possible
k = numel(files);
SylParam(k).name = '';
SylParam(k).z_Dur = [];
SylParam(k).z_FFreq1 = [];
SylParam(k).z_FFreq2 = [];
SylParam(k).z_FFreq3 = [];
SylParam(k).z_FFreq4 = [];
SylParam(k).z_HPAmp = [];
SylParam(k).z_FSlope1 = [];
SylParam(k).z_FSlope2 = [];
SylParam(k).z_FSlope3 = [];
SylParam(k).z_FSlope4 = [];
SylParam(k).z_ASlope12 = [];
SylParam(k).z_ASlope23 = [];
SylParam(k).z_ASlope34 = [];
SylParam(k).z_SEntropy1 = [];
SylParam(k).z_SEntropy2 = [];
SylParam(k).z_SEntropy3 = [];
SylParam(k).z_SEntropy4 = [];
SylParam(k).z_TEntropy1 = [];
SylParam(k).z_TEntropy2 = [];
SylParam(k).z_TEntropy3 = [];
SylParam(k).z_TEntropy4 = [];
SylParam(k).z_STEntropy1 = [];
SylParam(k).z_STEntropy2 = [];
SylParam(k).z_STEntropy3 = [];
SylParam(k).z_STEntropy4 = [];
SylParam(k).z_AM1 = [];
SylParam(k).z_AM2 = [];
SylParam(k).z_AM3 = [];
SylParam(k).z_AM4 = [];
SylParam(k).z_FM1 = [];
SylParam(k).z_FM2 = [];
SylParam(k).z_FM3 = [];
SylParam(k).z_FM4 = [];

raw_Params = zeros(numel(files),33);
z_Params = zeros(numel(files),33);

%Extract acoustic features from the syllable
for k=1:numel(files);
    SylParam(k).name = files(k).name;
    path = [dest '\' SylParam(k).name];
    Params = ParamExtract33(path)';  %[SylParam(k).Dur, SylParam(k).FFreq, SylParam(k).HPAmp, SylParam(k).FSlope, SylParam(k).ASlope, SylParam(k).SEntropy, SylParam(k).TEntropy, SylParam(k).STEntropy, SylParam(k).AM, SylParam(k).FM]=ParamExtract(path);
    raw_Params(k,1:33) = Params;
end

clear files;

%Calculate Z-scores for each of the spectral components captured
m_rawParams = mean(raw_Params);
std_rawParams = std(raw_Params);
[r c] = size(raw_Params);
for i=1:c
    z_Params(:,i) = (raw_Params(:,i)-m_rawParams(i))./std_rawParams(i);
end

%Copy all z-scores out to SylParam for permenant matching record
for k = 1:numel(SylParam)
    SylParam(k).z_Dur = z_Params(k,1);
    SylParam(k).z_FFreq1 = z_Params(k,2);
    SylParam(k).z_FFreq2 = z_Params(k,3);
    SylParam(k).z_FFreq3 = z_Params(k,4);
    SylParam(k).z_FFreq4 = z_Params(k,5);
    SylParam(k).z_HPAmp = z_Params(k,6);
    SylParam(k).z_FSlope1 = z_Params(k,7);
    SylParam(k).z_FSlope2 = z_Params(k,8);
    SylParam(k).z_FSlope3 = z_Params(k,9);
    SylParam(k).z_FSlope4 = z_Params(k,10);
    SylParam(k).z_ASlope12 = z_Params(k,11);
    SylParam(k).z_ASlope23 = z_Params(k,12);
    SylParam(k).z_ASlope34 = z_Params(k,13);
    SylParam(k).z_SEntropy1 = z_Params(k,14);
    SylParam(k).z_SEntropy2 = z_Params(k,15);
    SylParam(k).z_SEntropy3 = z_Params(k,16);
    SylParam(k).z_SEntropy4 = z_Params(k,17);
    SylParam(k).z_TEntropy1 = z_Params(k,18);
    SylParam(k).z_TEntropy2 = z_Params(k,19);
    SylParam(k).z_TEntropy3 = z_Params(k,20);
    SylParam(k).z_TEntropy4 = z_Params(k,21);
    SylParam(k).z_STEntropy1 = z_Params(k,22);
    SylParam(k).z_STEntropy2 = z_Params(k,23);
    SylParam(k).z_STEntropy3 = z_Params(k,24);
    SylParam(k).z_STEntropy4 = z_Params(k,25);
    SylParam(k).z_AM1 = z_Params(k,26);
    SylParam(k).z_AM2 = z_Params(k,27);
    SylParam(k).z_AM3 = z_Params(k,28);
    SylParam(k).z_AM4 = z_Params(k,29);
    SylParam(k).z_FM1 = z_Params(k,30);
    SylParam(k).z_FM2 = z_Params(k,31);
    SylParam(k).z_FM3 = z_Params(k,32);
    SylParam(k).z_FM4 = z_Params(k,33);
end

% [coeff score latent]=princomp(DataSet);
%  
% for i = 1:length(latent)
%     vartot(i)=sum(latent(1:i))/length(latent);
% end
% 
% num_pc = find(vartot>.95, 1); %Determine how many prin components are required for 95% variation

end %end function




