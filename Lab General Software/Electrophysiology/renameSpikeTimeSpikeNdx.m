function renameSpikeTimeSpikeNdx(num)
cell = loadCell(num);
for(nMotif = 1:length(cell.motifData))
    if(isfield(cell.motifData{nMotif}, 'spiketime'))
        cell.motifData{nMotif}.spikeNdx = cell.motifData{nMotif}.spiketime;
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'spiketime');
    end
    if(isfield(cell.motifData{nMotif}, 'spikeTime'))
        cell.motifData{nMotif}.spikeNdx = cell.motifData{nMotif}.spikeTime;
        cell.motifData{nMotif} = rmfield(cell.motifData{nMotif}, 'spikeTime');
    end
end
saveCell(cell);