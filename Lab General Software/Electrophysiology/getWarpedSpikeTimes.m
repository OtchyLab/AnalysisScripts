function spikeTime = getWarpedSpikeTimes(cell, nMotif, warpTemplate)

sampRate = cell.exper.desiredInSampRate;

motifData = cell.motifData{nMotif};
spikeNdx = motifData.spikeNdx - (motifData.marker - 1);
warpMarks = motifData.warpMarks;

spikeTime = ndx2time(spikeNdx, sampRate);
warpTime = ndx2time(warpMarks - (motifData.marker-1) , sampRate);
warpTemplate = ndx2time(warpTemplate, sampRate);

%figure;
%plot(warpTime, -ones(length(warpTime),1),'+');
%plot(spikeTime, -ones(length(spikeTime),1),'x');
ylim([0,40]);
hold on;
%Align start 
spikeTime = spikeTime - (warpTime(1) - warpTemplate(1));
warpTime = warpTime - (warpTime(1) - warpTemplate(1));


%plot(warpTemplate,zeros(length(warpTime),1),'o');
%plot(warpTime,ones(length(warpTime),1),'+'); 
%plot(spikeTime, ones(length(spikeTime),1),'x');
for(nWarp = 2:length(warpTime))
    scale = (warpTemplate(nWarp) - warpTemplate(nWarp-1)) / (warpTime(nWarp) - warpTime(nWarp-1));
    shift = (warpTemplate(nWarp) - warpTemplate(nWarp-1)) - (warpTime(nWarp) - warpTime(nWarp-1));
    
    %spikes in the current warp zone
    ndx = find(spikeTime > warpTime(nWarp-1) & spikeTime <= warpTime(nWarp));
    spikeTime(ndx) = ((spikeTime(ndx) - warpTime(nWarp-1))*scale) + warpTime(nWarp-1);
    warpTime(nWarp) =((warpTime(nWarp) - warpTime(nWarp-1))*scale) + warpTime(nWarp-1);
    
    %spikes not in the current warp zone
    ndx = find(spikeTime > warpTime(nWarp));
    spikeTime(ndx) = spikeTime(ndx) + shift;
    if(nWarp~=length(warpTime))
        warpTime(nWarp+1:end) = warpTime(nWarp+1:end) + shift; 
    end
    %plot(warpTime,nWarp*ones(length(warpTime),1),'+'); 
    %plot(spikeTime, nWarp*ones(length(spikeTime),1),'x');
end