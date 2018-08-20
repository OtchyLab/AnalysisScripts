function displayCellMotif(cell, motifNum, numXPlots, numYPlots, nAudioPlot, nSigPlot, bMarkSpikes)
%Last four arguments are options.
if(~exist('numXPlots'))
    numXPlots = 2;
end
if(~exist('numYPlots'))
    numYPlots = 1;
end
if(~exist('nAudioPlot'))
    nAudioPlot = 1;
end
if(~exist('nSigPlot'))
    nSigPlot = 2;
end
if(~exist('bMarkSpikes'))
    bMarkSpikes = true;
end

subplot(numXPlots,numYPlots,nAudioPlot);
displayMotifAudio(cell, motifNum);

subplot(numXPlots,numYPlots,nSigPlot);
displayMotifTrace(cell, motifNum,bMarkSpikes);
