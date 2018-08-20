function [IntMat] = IntervalMatrixDIR
%This script loads the dataset from the saved process file (from StretchEm)
%and returns:
%IntMat = m x n matrix of all interval lengths
%       m = renditions
%       n = interval (syllables and gaps)
%keys = the filesnames of each rendition
%startSyl = the index of the starting syllable of sequence for the file
%sequence = the syllable sequence represented in the matrix
%RemapIndx = the key for the dir/undir renditions
%Breaks = the transition index
%Differs from the previous version (IntervalMatrix) in that it parses the
%data for Dir/Undir

%Request and load dataset
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\DirUndir XCorr Sets\Alignment Sets\*.mat');
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

%Find the break point for the dir/und
RemapIndx = data.directRemapInx;
breaks = find(diff(RemapIndx)<0,1,'first')+1;

%Special case for this dataset!!!
%breaks = 31;

IntMat = IntMat(RemapIndx,:);

%Plot interval output stats
intervals.un_m = mean(IntMat(1:breaks-1,:),1);
intervals.un_std = std(IntMat(1:breaks-1,:));
intervals.un_cv = intervals.un_std./intervals.un_m;
intervals.un_cvS = mean(intervals.un_cv(1:2:end));
intervals.un_cvG = mean(intervals.un_cv(2:2:end));

intervals.dir_m = mean(IntMat(breaks:end,:),1);
intervals.dir_std = std(IntMat(breaks:end,:));
intervals.dir_cv = intervals.dir_std./intervals.dir_m;
intervals.dir_cvS = mean(intervals.dir_cv(1:2:end));
intervals.dir_cvG = mean(intervals.dir_cv(2:2:end));

fig1 = figure;
subplot(2,1,1)
hold on
bar((2:2:length(intervals.un_m)*2)-.25,intervals.un_m,'BarWidth',0.2)
bar((2:2:length(intervals.un_m)*2)+.25,intervals.dir_m,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar((2:2:length(intervals.un_m)*2)-.25,intervals.un_m,intervals.un_std,'.k');
errorbar((2:2:length(intervals.un_m)*2)+.25,intervals.dir_m,intervals.dir_std,'.k');
set(gca,'XTick',[],'XTickLabel',[])
set(gca,'TickDir','out','Box','off')
ylabel('Mean Int Duration (ms)')
xlim([1,length(intervals.un_m)*2+1])
tits = ['Interval lengths for ' fname(1:end-11)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
line([1,length(intervals.un_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar((2:2:length(intervals.un_m)*2),intervals.dir_m./intervals.un_m,'BarWidth',0.4,'FaceColor',[1 0 0])
errorbar((2:2:length(intervals.un_m)*2),intervals.un_m./intervals.un_m,intervals.un_std./intervals.un_m,'.k');
errorbar((2:2:length(intervals.un_m)*2),intervals.dir_m./intervals.un_m,intervals.dir_std./intervals.un_m,'.k');
temp = [];
nums = 1:.5:5;
for i = 1:2:length(intervals.un_cv)
    temp = [temp; char(num2str(nums(i)))];
    temp = [temp; char('G')];
end
set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
set(gca,'TickDir','out','Box','off')
ylabel('Norm Duration')
xlabel('Syllable/Gap')
xlim([1,length(intervals.un_m)*2+1])
ylim([0.9,1.1])
tits = ['Normalized interval lengths'];
title(tits,'Interpreter','none')

fig2 = figure;
subplot(2,1,1)
hold on
scatter((2:2:length(intervals.un_m)*2),intervals.un_cv,'sb')
scatter((2:2:length(intervals.un_m)*2),intervals.dir_cv,'sr')
set(gca,'TickDir','out','Box','off')
ylabel('CV')
xlim([1,length(intervals.un_m)*2+1])
tits = ['CV for ' fname(1:end-11)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
line([1,length(intervals.un_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar((2:2:length(intervals.un_m)*2),intervals.dir_cv./intervals.un_cv,'BarWidth',0.4,'FaceColor',[1 0 0])
set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
set(gca,'TickDir','out','Box','off')
xlabel('Syllable/Gap')
ylabel('Norm CV')
xlim([1,length(intervals.un_m)*2+1])
ylim([0.3,1.1])
title('Normalized CV')

% %Do the variability decomposition with Glaze and Troyer Code
% [intervals.un_model,intervals.un_z,intervals.un_u,intervals.un_eta,W0,sigma0,psi0] = timing_var_EM(IntMat(1:breaks-1,:),1);
% [intervals.dir_model,intervals.dir_z,intervals.dir_u,intervals.dir_eta,W0,sigma0,psi0] = timing_var_EM(IntMat(breaks:end,:),1);
% 
% fig3 = figure;
% subplot(2,3,1)
% hold on
% bar((2:2:length(intervals.un_m)*2)-.25,intervals.un_model.W,'BarWidth',0.2)
% bar((2:2:length(intervals.un_m)*2)+.25,intervals.dir_model.W,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.un_m)*2+1])
% title('Global Variability');
% 
% subplot(2,3,2)
% hold on
% bar((2:2:length(intervals.un_m)*2)-.25,intervals.un_model.psi,'BarWidth',0.2)
% bar((2:2:length(intervals.un_m)*2)+.25,intervals.dir_model.psi,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.un_m)*2+1])
% title('Independent Variability');
% 
% subplot(2,3,3)
% hold on
% bar((2:2:length(intervals.un_m)*2-2)-.25+1,intervals.un_model.sigma,'BarWidth',0.2)
% bar((2:2:length(intervals.un_m)*2-2)+.25+1,intervals.dir_model.sigma,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.un_m)*2+1])
% title('Jitter');
% 
% subplot(2,3,4:6)
% hold on
% bar((2:2:length(intervals.un_m)*2)-.25,intervals.dir_model.W./intervals.un_model.W,'BarWidth',0.2,'FaceColor',[1 0 0])
% bar((2:2:length(intervals.un_m)*2)+.25,intervals.dir_model.psi./intervals.un_model.psi,'BarWidth',0.2,'FaceColor',[0 1 0])
% bar((3:2:length(intervals.un_m)*2),intervals.dir_model.sigma./intervals.un_model.sigma,'BarWidth',0.2,'FaceColor',[0 0 1])
% line([1,length(intervals.un_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
% set(gca,'XTick',2:2:length(intervals.un_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlabel('Syllable/Gap')
% ylabel('Variance Normalized to Undirected')
% xlim([1,length(intervals.un_m)*2+1])
% %ylim([0.3,1.1])
% title('Deconstructed Variability (Global,Independent,Jitter)')

%Save data to file
savename = [pathLoc, 'Intervals\' fname(1:end-11), ' DIR ints.mat'];
save(savename,'IntMat','filenames','sequence','intervals','RemapIndx','breaks')

saveas(fig1,[pathLoc, 'Intervals\Fig\' fname(1:end-11), ' durations.fig'])
saveas(fig1,[pathLoc, 'Intervals\Tiff\' fname(1:end-11), ' durations.tif'])
saveas(fig2,[pathLoc, 'Intervals\Fig\' fname(1:end-11), ' CV.fig'])
saveas(fig2,[pathLoc, 'Intervals\Tiff\' fname(1:end-11), ' CV.tif'])

close(fig1)
close(fig2)