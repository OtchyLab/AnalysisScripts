function spikeSoundCorr(spikeCorr,simScore)

%input: spikeCorr: correlation matrix for spikes
%simScore: similarity score vector from SAP
simSize=sqrt(length(simScore));
if simSize~=size(spikeCorr,1)
    warndlg('two matrices not the same size')
    return;
end
for i=1:simSize
    simMatrix(i,1:simSize)=simScore(1+(i-1)*simSize:i*simSize);
end
spikeCorrVector=averageMatrix(spikeCorr,2);
soundCorrVector=averageMatrix(simMatrix,1);
syll_variability=mean(soundCorrVector)
figure;plot(spikeCorrVector,soundCorrVector,'.');
[R,P] = corrcoef(spikeCorrVector,soundCorrVector) 
