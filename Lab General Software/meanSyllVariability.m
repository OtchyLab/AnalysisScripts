function meanSyllVariability(simScore)

%input: spikeCorr: correlation matrix for spikes
%simScore: similarity score vector from SAP
simSize=sqrt(length(simScore));
for i=1:simSize
    simMatrix(i,1:simSize)=simScore(1+(i-1)*simSize:i*simSize);
end

soundCorrVector=averageMatrix(simMatrix,1);
syll_variability=mean(soundCorrVector)
