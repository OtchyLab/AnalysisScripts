function [IntMat] = IntervalMatrixStop
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
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\Stopping Analysis\*.mat');
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

%Create indices for each of the analysis groups
compMotif = find(data.postSyll==5);
stopMotif = find(data.postSyll==99);

%Plot interval output stats
intervals.comp_m = mean(IntMat(compMotif,:),1);
intervals.comp_std = std(IntMat(compMotif,:),0,1);
intervals.comp_cv = intervals.comp_std./intervals.comp_m;
intervals.comp_cvS = mean(intervals.comp_cv(1:2:end));
intervals.comp_cvG = mean(intervals.comp_cv(2:2:end));
intervals.comp_totalLen = sum(IntMat(compMotif,:),2);

intervals.stop_m = mean(IntMat(stopMotif,:),1);
intervals.stop_std = std(IntMat(stopMotif,:),0,1);
intervals.stop_cv = intervals.stop_std./intervals.stop_m;
intervals.stop_cvS = mean(intervals.stop_cv(1:2:end));
intervals.stop_cvG = mean(intervals.stop_cv(2:2:end));
intervals.stop_totalLen = sum(IntMat(stopMotif,:),2);

fig1 = figure;
subplot(2,1,1)
hold on
bar((2:2:length(intervals.comp_m)*2)-.25,intervals.comp_m,'BarWidth',0.2)
bar((2:2:length(intervals.comp_m)*2)+.25,intervals.stop_m,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar((2:2:length(intervals.comp_m)*2)-.25,intervals.comp_m,intervals.comp_std,'.k');
errorbar((2:2:length(intervals.comp_m)*2)+.25,intervals.stop_m,intervals.stop_std,'.k');
set(gca,'XTick',[],'XTickLabel',[])
set(gca,'TickDir','out','Box','off')
ylabel('Mean Int Duration (ms)')
xlim([1,length(intervals.comp_m)*2+1])
tits = ['Interval lengths for ' fname(1:end-11)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
line([1,length(intervals.comp_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar((2:2:length(intervals.comp_m)*2),intervals.stop_m./intervals.comp_m,'BarWidth',0.4,'FaceColor',[1 0 0])
errorbar((2:2:length(intervals.comp_m)*2),intervals.comp_m./intervals.comp_m,intervals.comp_std./intervals.comp_m,'.k');
errorbar((2:2:length(intervals.comp_m)*2),intervals.stop_m./intervals.comp_m,intervals.stop_std./intervals.comp_m,'.k');
temp = [];
nums = 1:.5:5;
for i = 1:2:length(intervals.comp_cv)
    temp = [temp; char(num2str(nums(i)))];
    temp = [temp; char('G')];
end
set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
set(gca,'TickDir','out','Box','off')
ylabel('Norm Duration')
xlim([1,length(intervals.comp_m)*2+1])
ylim([0.5,1.5])
tits = ['Normalized interval lengths'];
title(tits,'Interpreter','none')

fig2 = figure;
subplot(2,1,1)
hold on
scatter((2:2:length(intervals.comp_m)*2),intervals.comp_cv,'sb')
scatter((2:2:length(intervals.comp_m)*2),intervals.stop_cv,'sr')
set(gca,'TickDir','out','Box','off')
ylabel('CV')
xlim([1,length(intervals.comp_m)*2+1])
tits = ['CV for ' fname(1:end-11)];
title(tits,'Interpreter','none')

subplot(2,1,2)
hold on
line([1,length(intervals.comp_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar((2:2:length(intervals.comp_m)*2),intervals.stop_cv./intervals.comp_cv,'BarWidth',0.4,'FaceColor',[1 0 0])
set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
set(gca,'TickDir','out','Box','off')
xlabel('Syllable/Gap')
ylabel('Norm CV')
xlim([1,length(intervals.comp_m)*2+1])
ylim([0.25,1.25])
title('Normalized CV')

fig3 = figure;
subplot(2,1,1)
hold on
scatter(1*ones(length(intervals.comp_totalLen),1)+.1*randn(length(intervals.comp_totalLen),1),intervals.comp_totalLen,'.b')
scatter(3*ones(length(intervals.stop_totalLen),1)+.1*randn(length(intervals.stop_totalLen),1),intervals.stop_totalLen,'.r')
errorbar([1,3],[mean(intervals.comp_totalLen),mean(intervals.stop_totalLen)],[std(intervals.comp_totalLen),std(intervals.stop_totalLen)],'sk')
set(gca,'XTick',[1,3],'XTickLabel',{char('Complete'),char('Stopped')})
set(gca,'TickDir','out','Box','off')
ylabel('Motif Length (ms)')
xlim([0,4])
title('Motif Duration and Stopping')

subplot(2,1,2)
hold on
scatter(1*ones(length(intervals.comp_totalLen),1)+.1*randn(length(intervals.comp_totalLen),1),IntMat(compMotif,end),'xb')
scatter(3*ones(length(intervals.stop_totalLen),1)+.1*randn(length(intervals.stop_totalLen),1),IntMat(stopMotif,end),'xr')
errorbar([1,3],[mean(IntMat(compMotif,end)),mean(IntMat(stopMotif,end))],[std(IntMat(compMotif,end)),std(IntMat(stopMotif,end))],'sk')
set(gca,'XTick',[1,3],'XTickLabel',{char('Complete'),char('Stopped')})
set(gca,'TickDir','out','Box','off')
ylabel('Last Syllable Length (ms)')
xlim([0,4])
title('Pre-Target Syllable Duration and Stopping')

% %Do the variability decomposition with Glaze and Troyer Code
% [intervals.comp_model,intervals.comp_z,intervals.comp_u,intervals.comp_eta,W0,sigma0,psi0] = timing_var_EM(IntMat(1:breaks-1,:),1);
% [intervals.stop_model,intervals.stop_z,intervals.stop_u,intervals.stop_eta,W0,sigma0,psi0] = timing_var_EM(IntMat(breaks:end,:),1);
% 
% fig3 = figure;
% subplot(2,3,1)
% hold on
% bar((2:2:length(intervals.comp_m)*2)-.25,intervals.comp_model.W,'BarWidth',0.2)
% bar((2:2:length(intervals.comp_m)*2)+.25,intervals.stop_model.W,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.comp_m)*2+1])
% title('Global Variability');
% 
% subplot(2,3,2)
% hold on
% bar((2:2:length(intervals.comp_m)*2)-.25,intervals.comp_model.psi,'BarWidth',0.2)
% bar((2:2:length(intervals.comp_m)*2)+.25,intervals.stop_model.psi,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.comp_m)*2+1])
% title('Independent Variability');
% 
% subplot(2,3,3)
% hold on
% bar((2:2:length(intervals.comp_m)*2-2)-.25+1,intervals.comp_model.sigma,'BarWidth',0.2)
% bar((2:2:length(intervals.comp_m)*2-2)+.25+1,intervals.stop_model.sigma,'BarWidth',0.2,'FaceColor',[1 0 0])
% set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlim([1,length(intervals.comp_m)*2+1])
% title('Jitter');
% 
% subplot(2,3,4:6)
% hold on
% bar((2:2:length(intervals.comp_m)*2)-.25,intervals.stop_model.W./intervals.comp_model.W,'BarWidth',0.2,'FaceColor',[1 0 0])
% bar((2:2:length(intervals.comp_m)*2)+.25,intervals.stop_model.psi./intervals.comp_model.psi,'BarWidth',0.2,'FaceColor',[0 1 0])
% bar((3:2:length(intervals.comp_m)*2),intervals.stop_model.sigma./intervals.comp_model.sigma,'BarWidth',0.2,'FaceColor',[0 0 1])
% line([1,length(intervals.comp_m)*2+1],[1,1],'Color','b','LineStyle','--','LineWidth',1)
% set(gca,'XTick',2:2:length(intervals.comp_cv)*2,'XTickLabel',temp)
% set(gca,'TickDir','out','Box','off')
% xlabel('Syllable/Gap')
% ylabel('Variance Normalized to Undirected')
% xlim([1,length(intervals.comp_m)*2+1])
% %ylim([0.3,1.1])
% title('Deconstructed Variability (Global,Independent,Jitter)')

%Save data to file
savename = [pathLoc, 'Intervals\' fname(1:end-11), ' stopping ints.mat'];
save(savename,'IntMat','filenames','sequence','intervals','compMotif','stopMotif')

saveas(fig1,[pathLoc, 'Intervals\Fig\' fname(1:end-11), ' durations.fig'])
saveas(fig1,[pathLoc, 'Intervals\Tiff\' fname(1:end-11), ' durations.tif'])
saveas(fig2,[pathLoc, 'Intervals\Fig\' fname(1:end-11), ' CV.fig'])
saveas(fig2,[pathLoc, 'Intervals\Tiff\' fname(1:end-11), ' CV.tif'])
saveas(fig3,[pathLoc, 'Intervals\Fig\' fname(1:end-11), ' motif and pretarget.fig'])
saveas(fig3,[pathLoc, 'Intervals\Tiff\' fname(1:end-11), ' motif and pretarget.tif'])

close(fig1)
close(fig2)
close(fig3)
