function showCellFile(cellFile)

position=[10,50,920,1000];
figure('position',position,'PaperOrientation','landscape','PaperPositionMode','auto');
subplot('position',[0.07 0.7 0.3 0.2] );
lags=-1:0.001:1;
smoothedAC=smooth(cellFile.ACF,10);
plot(lags,smoothedAC);
xlabel('Time lag (s)');
ylabel('Correlation');
title('Auto correlation');
subplot('position',[0.07 0.1 0.3 0.2] );
% plot(cellFile.ISIDrug{1,2},cellFile.ISIDrug{1,1},'color','r');
% hold on;
% yy = smooth(cellFile.ISINoDrug{1,1},0.1,'rlowess');
% plot(cellFile.ISINoDrug{1,2},cellFile.ISINoDrug{1,1}, cellFile.ISINoDrug{1,2},yy);
% xlabel('Hz');
% ylabel('Relative Frequency');
% title('IFR Distribution');
hold off;
subplot('position',[0.07 0.4 0.3 0.2] );
 semilogy(cellFile.activeDrug{1,2},cellFile.activeDrug{1,1},'color','r');
%plot(cellFile.activeDrug{1,2},cellFile.activeDrug{1,1},'color','r');
hold on;
 semilogy(cellFile.activeNoDrug{1,2},cellFile.activeNoDrug{1,1});
plot(cellFile.activeNoDrug{1,2},cellFile.activeNoDrug{1,1});
xlabel('Hz');
ylabel('Fraction of time active');
title('Fraction active ');
hold off;

subplot('position',[0.5 0.4 0.3 0.2] );
for i=1:50 bins(i)=(0.0002)*1.15^i;
end
ISIHist=histc(cellFile.spikeI,bins)/length(cellFile.spikeI);
ISIHistDrugs=histc(cellFile.spikeIDrugs,bins)/length(cellFile.spikeIDrugs);
ISIHistCum = cumsum(ISIHist);
if (sum(ISIHistDrugs)>0)
    ISIHistDrugsCum=cumsum(ISIHistDrugs);
    loglog(bins,ISIHistDrugsCum,'r');
    hold on;
end
loglog(bins,ISIHistCum);
xlabel('ISI');
ylabel('Relative Frequency');
title('ISI distribution');
xlim([0.001 1]);
hold off;

% subplot('position',[0.5 0.55 0.3 0.35] );
% cellFile.frDrug{1,1}=cellFile.frDrug{1,1}/sum(cellFile.frDrug{1,1});
% if (isfield(cellFile,'frNoDrug'))
%     cellFile.frNoDrug{1,1}=cellFile.frNoDrug{1,1}/sum(cellFile.frNoDrug{1,1});
%     plot(cellFile.frNoDrug{1,2},cellFile.frNoDrug{1,1});
%     hold on
% end
% plot(cellFile.frDrug{1,2},cellFile.frDrug{1,1},'color','r');
% xlabel('Hz');
% ylabel('Relative Frequency');
% title('Firing rate Singing');
% hold off;

subplot('position',[0.5 0.1 0.3 0.2] );
plot(cellFile.spikeI(1:end-1),cellFile.spikeI(2:end),'.','markersize',0.05);
hold on;
plot(cellFile.spikeIDrugs(1:end-1),cellFile.spikeIDrugs(2:end),'.','markersize',0.05,'color','r');

xlim([0 0.3]);
ylim([0 0.3]);
xlabel('ISI(i)');
% line([0 0.01],[0.3 0.01],'color','k');
% line([0.01 0],[0.01 0.3],'color','k');
ylabel('ISI(i+1)');
title('ISI(i) vs ISI (i+1)');
hold off;

if (~isempty(cellFile.spontSpikesNoDrug))
    xhist=(0:1:100);
    frSpont=1./cellFile.spontSpikesNoDrug;
    histSpont=histc(frSpont,xhist);
    histSpontNoDrug=histSpont/sum(histSpont);
%     plot(xhist,histSpontNoDrug);
%     hold on;
end
if (~isempty(cellFile.spontSpikesDrug))
    xhist=(0:1:100);
    frSpont=1./cellFile.spontSpikesDrug;
    histSpont=histc(frSpont,xhist);
    histSpontDrug=histSpont/sum(histSpont);
%     plot(xhist,histSpontDrug, 'color','r');
%     hold on;
end
% xlabel('Hz');
% ylabel('Relative Frequency');
% title('Firing rate Spontaneous');
% hold off;


burst=0;
for i=1:length(cellFile.spikeI)-1
    if (cellFile.spikeI(i)<0.01 && cellFile.spikeI(i+1)<(0.01))
        burst=burst+1;
    end
end
fraction=burst/length(cellFile.spikeI)

burst=0;
for i=1:length(cellFile.spikeIDrugs)-1
    if (cellFile.spikeIDrugs(i)<0.01 && cellFile.spikeIDrugs(i+1)<(0.01))
        burst=burst+1;
    end
end
if length(cellFile.spikeIDrugs)>0
    fractionDrugs=burst/length(cellFile.spikeIDrugs);
else
    fractionDrugs=0;
end
subplot('position',[0.82 0.1 0.1 0.75]);
text(0.1,1,['Cell Number ' num2str(cellFile.cellNo)]);
text(0.1,0.92,[num2str(cellFile.dph), ' days old']);
text(0.1,0.84,[num2str(cellFile.depth), ' microns deep']);
text(0.1,0.76,['Electrode number: ' num2str(cellFile.channel)]);
text(0.1,0.68,['Fraction ND: ' num2str(fraction)]);
text(0.1,0.60,['Fraction D: ',num2str(fractionDrugs)]);
if (~isempty(cellFile.spontSpikesDrug))
    spontDrug=sum(xhist.*histSpontDrug);
    text(0.1,0.52,['Spont. FR (D): ',num2str(spontDrug, '%6.1f')]);
end
if (~isempty(cellFile.frDrug))
    frDrug=cellFile.frDrug;
    text(0.1,0.36,['FR (D): ',num2str(frDrug,'%6.1f')]);   
end
if (~isempty(cellFile.spontSpikesNoDrug))
    spontNoDrug=sum(xhist.*histSpontNoDrug);
    text(0.1,0.44,['Sp.FR (ND): ',num2str(spontNoDrug ,'%6.1f')]);
end
if (~isempty(cellFile.frNoDrug))
    frNoDrug=cellFile.frNoDrug;
    text(0.1,0.28,['FR (D): ',num2str(frNoDrug ,'%6.1f')]);
end
text(0.1,0.2,['Sequence ',num2str(cellFile.Seq)]);
text(0.1,0.12,['P-value for FR: ',num2str(cellFile.FRsignificance ,'%6.2f')]);
text(0.1,0.04,['PSTH Corr: ',num2str(cellFile.PSTHcorr)]);

axis off;

