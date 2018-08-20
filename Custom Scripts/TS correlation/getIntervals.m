function intervals = getIntervals()

 %load('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur639\Pur639 120326\Pur639_120326_dataset (platonic,Ch2).mat')

 %Select and load file containing dataset
[fname, loc] = uigetfile('*.mat', 'Select file containing dataset.');
if isempty(fname) || isempty(loc)
    set(handles.text_message,'String','No Dataset selected. Process cancelled.')
    return
end
load([loc fname]);
 
 av_spectrogram_day1=squeeze(mean(data.aligned_audioCube,1)); 
 av_neuroTS_day1=squeeze(mean(data.aligned_neuroTS.^2,1));
 
 intervals_day1 = [];
 a = transpose(data.templatesyllBreaks);
 
 for i = 1:83
    [warpedOut] = getWarpedStarts([data.p{i},data.q{i}],a(:));
    intervals_day1 = [intervals_day1;diff(warpedOut)'];
 end
 
 intervals = mean(intervals_day1)