function StitchTogether
%The script stitches together the audio and neural power traces that
%correspend to the completed and stopped motifs.  It expects (as selected
%input) the the output of the MUA_Stopping script.

loc = 'C:\Users\Tim\Desktop\Stopping Analysis\Aligned\';
cd(loc)
files = dir('*Ch1*.mat');
numFiles = length(files);

maxLen = 0;
pntr = 0;
for i = 1:numFiles
    load(files(i).name,'xs','audioMat','neuroMat','compMotif','stopMotif')
    audioM_comp{i} = nanmean(audioMat(compMotif,:),1);
    audioM_stop{i} = nanmean(audioMat(stopMotif,:),1);
    neuroM_comp{i} = nanmean(neuroMat(compMotif,:),1);
    neuroM_stop{i} = nanmean(neuroMat(stopMotif,:),1);
    xses{i} = xs;
    starts(i) = xs(1);
    maxLen = max(maxLen,length(xs));
end

%Create blank/NaN'd sets
audioMesh_c = NaN(numFiles,maxLen);
audioMesh_s = NaN(numFiles,maxLen);
neuroMesh_c = NaN(numFiles,maxLen);
neuroMesh_s = NaN(numFiles,maxLen);

%And layer everything into the NaN'd sets
longPre = min(starts);
for i = 1:numFiles
    ind = (starts(i)-longPre+1):(xses{i}(end)-longPre+1);
    audioMesh_c(i,ind) = audioM_comp{i};
    audioMesh_s(i,ind) = audioM_stop{i};
    neuroMesh_c(i,ind) = neuroM_comp{i};
    neuroMesh_s(i,ind) = neuroM_stop{i};
    z = find(ind==0,1);
end

xmax = 140;
ticks = [19,39,59,79,99,119,139];
tLabels = {char('-60');char('-40');char('-20');char('0');char('20');char('40');char('60')};

% figure
% mesh(audioMesh_c);
% hold on
% mesh(audioMesh_s);
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')
% 
% figure
% mesh(neuroMesh_c);
% hold on
% mesh(neuroMesh_s);
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')
% 
% figure
% subplot(2,1,1)
% imagesc(audioMesh_c);
% line([-longPre,-longPre],[1,numFiles],'Color','w','LineStyle','--','LineWidth',2)
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')
% hold on
% subplot(2,1,2)
% imagesc(audioMesh_s);
% line([-longPre,-longPre],[1,numFiles],'Color','w','LineStyle','--','LineWidth',2)
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')
% 
% figure
% subplot(2,1,1)
% imagesc(neuroMesh_c);
% line([-longPre,-longPre],[1,numFiles],'Color','w','LineStyle','--','LineWidth',2)
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')
% hold on
% subplot(2,1,2)
% imagesc(neuroMesh_s);
% line([-longPre,-longPre],[1,numFiles],'Color','w','LineStyle','--','LineWidth',2)
% xlim([0,xmax])
% set(gca,'XTick',ticks,'XTickLabel',tLabels)
% set(gca,'TickDir','out','Box','off')

figure
subplot(2,1,1)
hold on
for i = 1:size(audioMesh_c,1)
    plot(audioMesh_c(i,:),'b');
end
plot(mean(audioMesh_c,1),'k','LineWidth',2);
plot(mean(audioMesh_c,1)+std(audioMesh_c,0,1),':k','LineWidth',2);
plot(mean(audioMesh_c,1)-std(audioMesh_c,0,1),':k','LineWidth',2);
for i = 1:size(audioMesh_s,1)
    plot(audioMesh_s(i,:),'r');
end
plot(mean(audioMesh_s,1),'g','LineWidth',2);
plot(mean(audioMesh_s,1)+std(audioMesh_s,0,1),':g','LineWidth',2);
plot(mean(audioMesh_s,1)-std(audioMesh_s,0,1),':g','LineWidth',2);
line([79,79],[-80,-35],'Color','k','LineStyle','--','LineWidth',2)
ylim([-80,-35])
xlim([19,119])
set(gca,'XTick',[],'XTickLabel',[])
set(gca,'TickDir','out','Box','off')
ylabel('Power (dB)')
title('Song Power in Completed and Stopping')

subplot(2,1,2)
hold on
for i = 1:size(neuroMesh_c,1)
    plot(neuroMesh_c(i,:),'b');
end
plot(mean(neuroMesh_c,1),'k','LineWidth',2);
plot(mean(neuroMesh_c,1)+std(neuroMesh_c,0,1),':k','LineWidth',2);
plot(mean(neuroMesh_c,1)-std(neuroMesh_c,0,1),':k','LineWidth',2);
for i = 1:size(neuroMesh_s,1)
    plot(neuroMesh_s(i,:),'r');
end
plot(mean(neuroMesh_s,1),'g','LineWidth',2);
plot(mean(neuroMesh_s,1)+std(neuroMesh_s,0,1),':g','LineWidth',2);
plot(mean(neuroMesh_s,1)-std(neuroMesh_s,0,1),':g','LineWidth',2);
line([79,79],[0,.35],'Color','k','LineStyle','--','LineWidth',2)
xlim([19,119])
ylim([0,0.35])
set(gca,'XTick',ticks,'XTickLabel',tLabels)
set(gca,'TickDir','out','Box','off')
xlabel('Time (ms)')
ylabel('Power (V^2)')
title('HVC Power in Completed and Stopping')


