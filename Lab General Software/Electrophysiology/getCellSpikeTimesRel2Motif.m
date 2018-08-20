function spikeTimes = getCellSpikeTimesRel2Motif(cell, warpTemplate)

sampRate = cell.exper.desiredInSampRate;
bWarp = exist('warpTemplate') && length(warpTemplate)>0;
    
spikeTimes = [];
for(nMotif = 1:length(cell.motifData))
    %Convert raw spike sample numbers in motifData.spikeNdx to
    %time in seconds relative to marker.
    motifData = cell.motifData{nMotif};
    if(isfield(motifData, 'spikeNdx'))
        marker = motifData.marker;       
        if(~bWarp)
            spikeTimeRelMarker = (motifData.spikeNdx - motifData.marker) / sampRate;
        else
            spikeTimeRelMarker = getWarpedSpikeTimes(cell, nMotif, warpTemplate);
        end
        
        spikeTimes = [spikeTimes, spikeTimeRelMarker];        
    end
end

