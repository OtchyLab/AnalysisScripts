function showConsqTraces(cell, startMotif)
h = figure;
cell = loadCell(cell.num);
subplot(2,1,1);
curr = 0;
for(nMotif = startMotif:length(cell.motifData))
    motifData = cell.motifData{nMotif};
    if(isfield(motifData,'spikeNdx'))
        subplot(2,1,curr+1);
        displayMotifTrace(cell,nMotif, false);
        title(num2str(nMotif));
        curr = mod(curr+1,2);
        pause;
        figure(h);
    end
end
