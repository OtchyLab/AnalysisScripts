function spikeTime = spikeTimeFromNdx(cell,nMotif)
%Takes the spikeNdx of a motif and converts it to spike times in seconds
%relative to motif marker.

spikeTime = cell.motifData{nMotif}.spikeNdx;
spikeTime = spikeTime - cell.motifData{nMotif}.marker;
spikeTime = spikeTime / cell.exper.desiredInSampRate;
