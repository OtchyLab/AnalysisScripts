function [nBinCounts, binEdges] = PSTHacrossCells(cellNums, cellDir, startTime, endTime, binSize)

spikeTimes = [];
for(cellNum = cellNums)
    cell = loadCell(cellNum, cellDir);
    spikeTimes = [spikeTimes, getCellSpikeTimesRel2Motif(cell)];
end

[nBinCounts, binEdges] = computePSTH(spikeTimes, startTime, endTime, binSize);
bar(binEdges,nBinCounts,'histc')