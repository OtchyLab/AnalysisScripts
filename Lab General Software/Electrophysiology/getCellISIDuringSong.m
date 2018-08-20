function isi = getCellISIDuringSong(cell)

sampRate = cell.exper.desiredInSampRate;
isi = [];
for(nMotif = 1:length(cell.motifData))
    %Convert raw spike sample numbers in motifData.spikeNdx to
    %time in seconds relative to marker.
    motifData = cell.motifData{nMotif};
    if(isfield(motifData, 'spikeNdx'))
        spikeTimeRelMarker = (motifData.spikeNdx - motifData.marker) / sampRate;
        isi = [isi, diff(spikeTimeRelMarker)];
    end
end