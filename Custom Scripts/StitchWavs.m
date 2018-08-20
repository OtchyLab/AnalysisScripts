function [Stitched] = StitchWavs(folder, gap)
%Take input folder and concatenate all .wav files in that folder into a
%single .wav file.  gap designates the number of milliseconds to insert
%between snippets.

%Sampling rate for input and output wavs
fs = 44100;

if nargin == 1
    gap = 15;
end

if nargin == 0
    gap = 15;
    [folder]=uigetdir('Song Directory');
end

%Convert gap (in ms) to sample counts
gap_counts = floor((gap/1000)*fs);
gap_segment = zeros(gap_counts,1);

%Get the list of .wavs to concatenate
files = dir([folder '\*.wav']);

Stitched = gap_segment;

for i=1:length(files)
    work_wav = wavread([folder '\' files(i).name]);
    Stitched = [Stitched;work_wav;gap_segment];
end



end