function cellRaster(cell)

sampRate = cell.exper.desiredInSampRate;
motifData = cell.motifData{1};
timeBuff = .1;
startTime = - motifData.marker / sampRate;
endTime = startTime + (motifData.ndxStop - motifData.ndxStart)/sampRate; 

%Plot raster of spikes from all motifs
hold on;
count = 1;
for(nMotif = 1:length(cell.motifData))
    %Convert raw spike sample numbers in motifData.spikeNdx to
    %time in seconds relative to marker.
    motifData = cell.motifData{nMotif};
    if(isfield(motifData, 'spikeNdx'))
        spikeTimeRelMarker = (motifData.spikeNdx - motifData.marker) / sampRate;
        plot(spikeTimeRelMarker, repmat(count,1,length(spikeTimeRelMarker)),'.');
        count = count + 1;
    end
end

xlim([startTime-timeBuff, endTime + timeBuff]);
title(['cell:', num2str(cell.num), ' Raster'])
xlabel('time secs');
ylabel('motif number (not counting motifs without spikes)');