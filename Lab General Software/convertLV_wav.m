function convertLV_wav(directory)
%fs is the sampling rate.

%syllStartTimes: the start time of each syllable.  This is not the onset of
%sound, but the time at which the audio level crosses above threshold.
%syllEndTimes: the end time of each syllable.  This is not the offset of
%sound, but the time at which the audio level crosses below threshold.
%noiseEst, noiseStd: the estimated noise level and variance.
%soundEst, soundStd: the estimated signal level and variance.

cd(directory);
d = dir ('_d*');    
fs=44100;
mkdir('wav files');

bDebug =false;
if ~isempty(d)
    for i=1:length(d)
        filename = d(i).name;
        fid = fopen(filename,'r','s');
        header=fread(fid,3,'single');
        audio=fread(fid,inf,'double');
        path=[directory,'\song files\',filename]
        wavwrite(audio,fs,path);

    end
end
