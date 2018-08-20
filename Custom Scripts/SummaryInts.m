function SummaryInts
%The script extracts the interesting features from the interval variables

loc = 'C:\Users\Tim\Desktop\DirUndir XCorr Sets\Alignment Sets\Intervals\';
cd(loc)
files = dir('*ints.mat');
numFiles = length(files);

maxLen = 0;
pntr = 0;
b1u = [];
b2u = [];
b3u = [];
b4u = [];
b1d = [];
b2d = [];
b3d = [];
b4d = [];
Allb1u = [];
Allb2u = [];
Allb3u = [];
Allb4u = [];
Allb1d = [];
Allb2d = [];
Allb3d = [];
Allb4d = [];
for i = 1:numFiles
    %Load the structure from the file
    load(files(i).name,'intervals','IntMat','breaks')
    
    %Sort the data by bird
    if i == 1 || i == 2 || i == 3
        b1u = [b1u; intervals.un_m];
        b1d = [b1d; intervals.dir_m];
        Allb1u = [Allb1u; IntMat(1:breaks-1,:)];
        Allb1d = [Allb1d; IntMat(breaks:end,:)];
    elseif i == 4 || i == 5
        b2u = [b2u; intervals.un_m];
        b2d = [b2d; intervals.dir_m];
        Allb2u = [Allb2u; IntMat(1:breaks-1,:)];
        Allb2d = [Allb2d; IntMat(breaks:end,:)];
    elseif i == 6 || i == 7 || i == 8
        b3u = [b3u; intervals.un_m];
        b3d = [b3d; intervals.dir_m];
        Allb3u = [Allb3u; IntMat(1:breaks-1,:)];
        Allb3d = [Allb3d; IntMat(breaks:end,:)];
    elseif i == 9 || i == 10
        b4u = [b4u; intervals.un_m];
        b4d = [b4d; intervals.dir_m];
        Allb4u = [Allb4u; IntMat(1:breaks-1,:)];
        Allb4d = [Allb4d; IntMat(breaks:end,:)];
    end
end

%Global Stretch all motifs to the un or dir man length
motifLen = sum(Allb1u,2);
Allb1u = Allb1u.*((motifLen/mean(motifLen))*ones(1,size(Allb1u,2)));
motifLen = sum(Allb1d,2);
Allb1d = Allb1d.*((motifLen/mean(motifLen))*ones(1,size(Allb1d,2)));

motifLen = sum(Allb2u,2);
Allb2u = Allb2u.*((motifLen/mean(motifLen))*ones(1,size(Allb2u,2)));
motifLen = sum(Allb2d,2);
Allb2d = Allb2d.*((motifLen/mean(motifLen))*ones(1,size(Allb2d,2)));

motifLen = sum(Allb3u,2);
Allb3u = Allb3u.*((motifLen/mean(motifLen))*ones(1,size(Allb3u,2)));
motifLen = sum(Allb3d,2);
Allb3d = Allb3d.*((motifLen/mean(motifLen))*ones(1,size(Allb3d,2)));

motifLen = sum(Allb4u,2);
Allb4u = Allb4u.*((motifLen/mean(motifLen))*ones(1,size(Allb4u,2)));
motifLen = sum(Allb4d,2);
Allb4d = Allb4d.*((motifLen/mean(motifLen))*ones(1,size(Allb4d,2)));

% b1delta = mean(b1d,1)./mean(b1u,1);
% b2delta = mean(b2d,1)./mean(b2u,1);
% b3delta = mean(b3d,1)./mean(b3u,1);
% b4delta = mean(b4d,1)./mean(b4u,1);
% 
% b1ds = std(b1d,0,1)./mean(b1u,1);
% b2ds = std(b2d,0,1)./mean(b2u,1);
% b3ds = std(b3d,0,1)./mean(b3u,1);
% b4ds = std(b4d,0,1)./mean(b4u,1);
% 
% b1us = std(b1u,0,1)./mean(b1u,1);
% b2us = std(b2u,0,1)./mean(b2u,1);
% b3us = std(b3u,0,1)./mean(b3u,1);
% b4us = std(b4u,0,1)./mean(b4u,1);

b1delta = mean(Allb1d,1)./mean(Allb1u,1);
b2delta = mean(Allb2d,1)./mean(Allb2u,1);
b3delta = mean(Allb3d,1)./mean(Allb3u,1);
b4delta = mean(Allb4d,1)./mean(Allb4u,1);

b1ds = std(Allb1d,0,1)./mean(Allb1u,1);
b2ds = std(Allb2d,0,1)./mean(Allb2u,1);
b3ds = std(Allb3d,0,1)./mean(Allb3u,1);
b4ds = std(Allb4d,0,1)./mean(Allb4u,1);

b1us = std(Allb1u,0,1)./mean(Allb1u,1);
b2us = std(Allb2u,0,1)./mean(Allb2u,1);
b3us = std(Allb3u,0,1)./mean(Allb3u,1);
b4us = std(Allb4u,0,1)./mean(Allb4u,1);

b1dcv = std(Allb1d,0,1)./mean(Allb1d,1);
b2dcv = std(Allb2d,0,1)./mean(Allb2d,1);
b3dcv = std(Allb3d,0,1)./mean(Allb3d,1);
b4dcv = std(Allb4d,0,1)./mean(Allb4d,1);

b1ucv = std(Allb1u,0,1)./mean(Allb1u,1);
b2ucv = std(Allb2u,0,1)./mean(Allb2u,1);
b3ucv = std(Allb3u,0,1)./mean(Allb3u,1);
b4ucv = std(Allb4u,0,1)./mean(Allb4u,1);

% Syl_delta = mean([b1delta(1:2:end),b2delta(1:2:end),b3delta(1:2:end),b4delta(1:2:end)]);
% Gap_delta = mean([b1delta(2:2:end),b2delta(2:2:end),b3delta(2:2:end),b4delta(2:2:end)]);
% Syl_s = std([b1delta(1:2:end),b2delta(1:2:end),b3delta(1:2:end),b4delta(1:2:end)]);
% Gap_s = std([b1delta(2:2:end),b2delta(2:2:end),b3delta(2:2:end),b4delta(2:2:end)]);

Syl_ucv = [b1ucv(1:2:end),b2ucv(1:2:end),b3ucv(1:2:end),b4ucv(1:2:end)];
Gap_ucv = [b1ucv(2:2:end),b2ucv(2:2:end),b3ucv(2:2:end),b4ucv(2:2:end)];
Syl_dcv = [b1dcv(1:2:end),b2dcv(1:2:end),b3dcv(1:2:end),b4dcv(1:2:end)];
Gap_dcv = [b1dcv(2:2:end),b2dcv(2:2:end),b3dcv(2:2:end),b4dcv(2:2:end)];

Syl_perCV = (Syl_dcv./Syl_ucv)-1;
Gap_perCV = (Gap_dcv./Gap_ucv)-1;

temp = [];
nums = 1:.5:5;
for i = 1:2:10
    temp = [temp; char(num2str(nums(i)))];
    temp = [temp; char('G')];
end

%Plot the output
figure
subplot(2,4,1)
hold on
line([1,length(b1delta)],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b1delta),b1delta,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar(1:length(b1delta),ones(length(b1delta),1),b1us,'.k')
errorbar(1:length(b1delta),b1delta,b1ds,'.k')
set(gca,'XTick',1:length(b1delta),'XTickLabel',temp(1:length(b1delta)))
set(gca,'TickDir','out','Box','off')
ylabel('Norm Duration')
xlim([0,length(b1delta)+1])
ylim([0.85,1.1])

subplot(2,4,2)
hold on
line([1,length(b2delta)],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b2delta),b2delta,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar(1:length(b2delta),ones(length(b2delta),1),b2us,'.k')
errorbar(1:length(b2delta),b2delta,b2ds,'.k')
set(gca,'XTick',1:length(b2delta),'XTickLabel',temp(1:length(b2delta)))
set(gca,'TickDir','out','Box','off')
%ylabel('Norm Duration')
%xlabel('Syllable/Gap')
xlim([0,length(b2delta)+1])
ylim([0.85,1.1])

subplot(2,4,3)
hold on
line([1,length(b3delta)],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b3delta),b3delta,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar(1:length(b3delta),ones(length(b3delta),1),b3us,'.k')
errorbar(1:length(b3delta),b3delta,b3ds,'.k')
set(gca,'XTick',1:length(b3delta),'XTickLabel',temp(1:length(b3delta)))
set(gca,'TickDir','out','Box','off')
xlim([0,length(b3delta)+1])
ylim([0.85,1.1])

subplot(2,4,4)
hold on
line([1,length(b4delta)],[1,1],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b4delta),b4delta,'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar(1:length(b4delta),ones(length(b4delta),1),b4us,'.k')
errorbar(1:length(b4delta),b4delta,b4ds,'.k')
set(gca,'XTick',1:length(b4delta),'XTickLabel',temp(1:length(b4delta)))
set(gca,'TickDir','out','Box','off')
xlim([0,length(b1delta)+1])
ylim([0.85,1.1])

ylims = [-.5,0.5];
subplot(2,4,5)
hold on
line([1,length(b1ucv)],[0,0],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b1ucv),(b1dcv./b1ucv)-1,'BarWidth',0.2,'FaceColor',[1 0 0])
set(gca,'XTick',1:length(b1delta),'XTickLabel',temp(1:length(b1delta)))
set(gca,'TickDir','out','Box','off')
ylabel('% change in CV from UNDIR')
xlim([0,length(b1ucv)+1])
ylim(ylims)

subplot(2,4,6)
hold on
line([1,length(b2ucv)],[0,0],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b2ucv),(b2dcv./b2ucv)-1,'BarWidth',0.2,'FaceColor',[1 0 0])
set(gca,'XTick',1:length(b2delta),'XTickLabel',temp(1:length(b2delta)))
set(gca,'TickDir','out','Box','off')
xlim([0,length(b2ucv)+1])
ylim(ylims)

subplot(2,4,7)
hold on
line([1,length(b3ucv)],[0,0],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b3ucv),(b3dcv./b3ucv)-1,'BarWidth',0.2,'FaceColor',[1 0 0])
set(gca,'XTick',1:length(b3delta),'XTickLabel',temp(1:length(b3delta)))
set(gca,'TickDir','out','Box','off')
xlim([0,length(b3ucv)+1])
ylim(ylims)

subplot(2,4,8)
hold on
line([1,length(b4ucv)],[0,0],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b4ucv),(b4dcv./b4ucv)-1,'BarWidth',0.2,'FaceColor',[1 0 0])
set(gca,'XTick',1:length(b4delta),'XTickLabel',temp(1:length(b4delta)))
set(gca,'TickDir','out','Box','off')
xlim([0,length(b4ucv)+1])
ylim(ylims)

figure
subplot(1,2,2)
hold on
line([1,3],[0,0],'Color','b','LineStyle','--','LineWidth',1)
scatter(ones(length(Syl_perCV),1)+.1*randn(length(Syl_perCV),1),Syl_perCV,'or')
scatter(3*ones(length(Gap_perCV),1)+.1*randn(length(Gap_perCV),1),Gap_perCV,'or')
%bar([1,3],[Syl_delta,Gap_delta],'BarWidth',0.2,'FaceColor',[1 0 0])
errorbar([1,3],[mean(Syl_perCV),mean(Gap_perCV)],[std(Syl_perCV),std(Gap_perCV)],'sk')
set(gca,'XTick',[1,3],'XTickLabel',{char('Syllables'),char('Gaps')})
set(gca,'YTick',[-0.5,0,0.5],'YTickLabel',{char('-0.5'),char('0'),char('0.5')})
set(gca,'TickDir','out','Box','off')
xlim([0,4])
ylim(ylims)
ylabel('Norm CV')

subplot(1,2,1)
hold on
scatter(Syl_dcv,Syl_ucv,'xb')
scatter(Gap_dcv,Gap_ucv,'sb')
line([0,.3],[0,.3],'Color','k','LineStyle','--','LineWidth',1)
set(gca,'XTick',[0,0.3],'XTickLabel',{char('0'),char('0.3')})
set(gca,'YTick',[0,0.3],'YTickLabel',{char('0'),char('0.3')})
set(gca,'TickDir','out','Box','off')
xlim([0,.3])
ylim([0,.3])
legend({'Syllables';'Gaps'})
xlabel('Directed CV')
ylabel('Undirected CV')








