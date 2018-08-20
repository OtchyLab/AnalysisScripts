function SummaryIntsStop
%The script extracts the interesting features from the interval variables

loc = 'C:\Users\Tim\Desktop\Stopping Analysis\Intervals\';
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

Allb1d = [];

for i = 1:numFiles
    %Load the structure from the file
    load(files(i).name,'intervals','IntMat','compMotif','stopMotif')

        b1u = [b1u; intervals.comp_m];
        b1d = [b1d; intervals.stop_m];
        Allb1u = [Allb1u; IntMat(compMotif,:)];
        Allb1d = [Allb1d; IntMat(stopMotif,:)];
end

%Global Stretch all motifs to the un or dir man length
motifLen = sum(Allb1u,2);
Allb1u = Allb1u.*((motifLen/mean(motifLen))*ones(1,size(Allb1u,2)));
motifLen = sum(Allb1d,2);
Allb1d = Allb1d.*((motifLen/mean(motifLen))*ones(1,size(Allb1d,2)));



b1delta = mean(Allb1d,1)./mean(Allb1u,1);

b1ds = std(Allb1d,0,1)./mean(Allb1u,1);


b1us = std(Allb1u,0,1)./mean(Allb1u,1);


b1dcv = std(Allb1d,0,1)./mean(Allb1d,1);


b1ucv = std(Allb1u,0,1)./mean(Allb1u,1);


Syl_ucv = [b1ucv(1:2:end)];
Gap_ucv = [b1ucv(2:2:end)];
Syl_dcv = [b1dcv(1:2:end)];
Gap_dcv = [b1dcv(2:2:end)];

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
subplot(2,1,1)
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

ylims = [-.5,0.5];
subplot(2,1,2)
hold on
line([1,length(b1ucv)],[0,0],'Color','b','LineStyle','--','LineWidth',1)
bar(1:length(b1ucv),(b1dcv./b1ucv)-1,'BarWidth',0.2,'FaceColor',[1 0 0])
set(gca,'XTick',1:length(b1delta),'XTickLabel',temp(1:length(b1delta)))
set(gca,'TickDir','out','Box','off')
ylabel('% change in CV from Comp Motif')
xlim([0,length(b1ucv)+1])
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
xlabel('Truncated CV')
ylabel('Completed CV')








