function cell = deleteSpikeTimes(cell)
cell = loadCell(cell.num);
for(nMotif = 1:length(cell.motifData))
    if(isfield(cell.motifData{nMotif}, 'spikeNdx'))
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'spikeNdx');
    end
    if(isfield(cell.motifData{nMotif}, 'bSpikesHandChecked'))
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'bSpikesHandChecked');
    end
    if(isfield(cell.motifData{nMotif}, 'spikeMethod'))
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'spikeMethod');
    end
    if(isfield(cell.motifData{nMotif}, 'spikeThreshold'))
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'spikeThreshold');
    end
end
saveCell(cell);