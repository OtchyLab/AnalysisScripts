function edges = getFreshEdges(Spec,starts)

%Collapse into power timeseries
specPowTS = -10*log10(sum(Spec));

%Use EM to fit 2 gaussians to distribution
gMix = gmdistribution.fit(specPowTS',2);
[powerM,pnt] = sort(gMix.mu);
temp = sqrt(squeeze(gMix.Sigma));
powerStd = temp(pnt);
thresh = powerM(1)+4*powerStd(1); %2SDs above the silent Gaussian mean

%Find the crossing points to estimate the onset and offsets of
%motifs
ind = 2:length(specPowTS);
Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);

onStarts = starts(1:2:end);
offStarts = starts(2:2:end);

edges = [];
 for k = 1:length(onStarts)
     [offset, pntr] = min(abs(Crossings.Up-onStarts(k)));
     edges = [edges;Crossings.Up(pntr)];
     
     [offset, pntr] = min(abs(Crossings.Down-offStarts(k)));
     edges = [edges;Crossings.Down(pntr)];
 end

 
figure
imagesc(-1*Spec)
hold on
for i = 1:length(edges)
    line([edges(i), edges(i)],[0,size(Spec,1)],'Color','k','LineWidth',2)
end
axis tight; axis xy;
title('Derived Template for Alignment')
xlabel('Time (ms)')
ylabel('Freq Bins')
 