function [nBinCounts, binEdges] = computePSTH(spikeTimes, startTime, endTime, binSize)

binEdges = startTime:binSize:endTime;
nBinCounts = histc(spikeTimes, binEdges)

%bar(binEdges,nBinCounts,'histc')

