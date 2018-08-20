function isi = getCellISINotSinging(cell)

sampRate = cell.exper.desiredInSampRate;
isi = [];
for(nNonSong = 1:length(cell.nonSongData))
    %Convert raw spike sample numbers in motifData.spikeNdx to
    %time in seconds relative to marker.
    nonSongData = cell.nonSongData{nNonSong};
    if(isfield(nonSongData, 'spikeNdx'))
        spikeTimeRelMarker = nonSongData.spikeNdx / sampRate;
        isi = [isi, diff(spikeTimeRelMarker)];
    end
end