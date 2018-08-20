function spikeTimes = warpSpikeTimes(spikeTimes, markTimes, templateMarkTimes)
%Aaron Andalman 2005
%Inefficient code for warping spikes, but plenty fast really.
%spikeTimes:  the times of the spikes relative to an origin within the
%motif.
%markTimes:  the warping points within this spikeTrain.
%templateMarkTimes:  the cononical times for the warping points.


%figure;
%plot(markTimes, -ones(length(markTimes),1),'+');
%plot(spikeTimes, -ones(length(spikeTimes),1),'x');
%ylim([0,40]);
%hold on;

%Align the starts
spikeTimes = spikeTimes - (markTimes(1) - templateMarkTimes(1));
markTimes = markTimes - (markTimes(1) - templateMarkTimes(1));

%plot(templateMarkTimes,zeros(length(markTimes),1),'o');
%plot(markTimes,ones(length(markTimes),1),'+'); 
%plot(spikeTimes, ones(length(spikeTimes),1),'x');

%Loop through each linear segment and warp spikes within this segment, and
%shift spikes after this segment.
for(nWarp = 2:length(markTimes))
    scale = (templateMarkTimes(nWarp) - templateMarkTimes(nWarp-1)) / (markTimes(nWarp) - markTimes(nWarp-1));
    shift = (templateMarkTimes(nWarp) - templateMarkTimes(nWarp-1)) - (markTimes(nWarp) - markTimes(nWarp-1));
    
    %spikes in the current warp zone
    ndx = find(spikeTimes > markTimes(nWarp-1) & spikeTimes <= markTimes(nWarp));
    spikeTimes(ndx) = ((spikeTimes(ndx) - markTimes(nWarp-1))*scale) + markTimes(nWarp-1);
    markTimes(nWarp) =((markTimes(nWarp) - markTimes(nWarp-1))*scale) + markTimes(nWarp-1);
    
    %spikes not in the current warp zone
    ndx = find(spikeTimes > markTimes(nWarp));
    spikeTimes(ndx) = spikeTimes(ndx) + shift;
    if(nWarp~=length(markTimes))
        markTimes(nWarp+1:end) = markTimes(nWarp+1:end) + shift; 
    end
    
    %plot(markTimes,nWarp*ones(length(markTimes),1),'+'); 
    %plot(spikeTimes, nWarp*ones(length(spikeTimes),1),'x');
end