function cell = addTimeWarpMarkers(cell, motifNums)
%motifNums is optional.  It can include the indices of the motifs you would
%like to add or replace time warp information for.

cell = loadCell(cell.num);
sampRate = cell.exper.desiredInSampRate;

if(~exist('motifNums'))
    motifNums = [1:length(cell.motifData)];
end

for(nMotif = motifNums)
    motif = cell.motifData{nMotif}
    marker = motif.marker;
    markerTime = marker/sampRate;
    
    %compute log power of audio
    logpow = log(motif.audio.^2);
    
    %Display audio
    time = [1:length(audio)] /sampRate;
    figure(1001);
    subplot(3,1,1);
    displayAudioSpecgram(audio, sampRate); axis tight;
    h1(1) = line([markerTime,markerTime],ylim); set(h1(1),'Color','red');
    subplot(3,1,2);
    plot(time, audio); axis tight;
    h2(1) = line([markerTime,markerTime],ylim); set(h2(1),'Color','red');
    subplot(3,1,3);
    plot(time, logpow); axis tight;
    h3(1) = line([markerTime,markerTime],ylim); set(h3(1),'Color','red');
    
    while(true)        
        
        
        [approxWarpMark, threshold, opt] = ginput(1)
        
        if(opt = 