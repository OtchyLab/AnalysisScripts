function [IntMat] = IntervalMatrix
%This script loads the dataset from the saved process file (from StretchEm)
%and returns:
%IntMat = m x n matrix of all interval lengths
%       m = renditions
%       n = interval (syllables and gaps)
%keys = the filesnames of each rendition
%startSyl = the index of the starting syllable of sequence for the file
%sequence = the syllable sequence represented in the matrix

%Request and load dataset
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\*.mat');
tots = [pathLoc, fname];
load(tots,'data','filenames','sequence')

%Select target syllable bounds from day 1 data
a = transpose(data.templatesyllBreaks);
a = a(:);

%Interval Durations from the paths and template files
rendNum = length(data.audio);
IntMat = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path = [data.p{i},data.q{i}];
    rendLength = size(data.audioSpecs_LN{i},2);
    
    %Calculate intervals and add to the stack
    [warpedOut] = getWarpedStarts(path,a(:));
    IntMat = [IntMat; diff(warpedOut)'];
    
end

%Plot interval output stats
intervals.m = mean(IntMat,1);
intervals.std = std(IntMat);
intervals.cv = intervals.std./intervals.m;

fig1 = figure;
subplot(2,1,1)
hold on
scatter(1:length(intervals.m),intervals.m,'or');
errorbar(1:length(intervals.m),intervals.m,intervals.std,'r');
set(gca,'XTick',[],'XTickLabel',[])
ylabel('Mean Int Duration (ms)')
xlim([0,length(intervals.m)+1])
tits = ['Interval lengths for ' fname(1:17)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
scatter(1:length(intervals.cv),intervals.cv,'sk')
line([0,length(intervals.cv)+1],[mean(intervals.cv(1:2:length(intervals.cv))),mean(intervals.cv(1:2:length(intervals.cv)))],'Color','r','LineStyle','--')
line([0,length(intervals.cv)+1],[mean(intervals.cv(2:2:length(intervals.cv))),mean(intervals.cv(2:2:length(intervals.cv)))],'Color','k','LineStyle','--')
temp = [];
nums = 1:.5:5;
for i = 1:2:length(intervals.cv)
    temp = [temp; char(num2str(nums(i)))];
    temp = [temp; char('G')];
end
set(gca,'XTick',1:length(intervals.cv),'XTickLabel',temp)
xlabel('Syllable/Gap')
ylabel('CV')
xlim([0,length(intervals.m)+1])

%Save data to file
% savename = [pathLoc, 'Intervals\' fname(1:17), ' ints.mat'];
% save(savename,'IntMat','filenames','sequence','intervals')
% 
% saveas(fig1,[pathLoc, 'Intervals\' fname(1:17), '.fig'])
%saveas(fig1,[pathLoc, 'Intervals\' fname(1:end-11), '.eps'])