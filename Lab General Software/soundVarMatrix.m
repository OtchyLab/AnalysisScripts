%function spikeSoundCorr(spikeCorr,simScore)
function soundVarMatrix(simScore)


%input: spikeCorr: correlation matrix for spikes
%simScore: similarity score vector from SAP
simSize=sqrt(length(simScore));
% if simSize~=size(spikeCorr,1)
%     warndlg('two matrices not the same size')
%     return;
% end
simScore=1-(100-simScore)/50;
for i=1:simSize
    simMatrix(i,1:simSize)=simScore(1+(i-1)*simSize:i*simSize);
end
figure;
simScore=sort(simScore);
min_sim=simScore(find(simScore>0.03,1,'first'));
max_sim=simScore(find(simScore<0.95,1,'last'));
imagesc([1 simSize],[simSize 1],simMatrix,[min_sim-0.05 max_sim]);
% spikeCorrVector=averageMatrix(spikeCorr,2);
% soundCorrVector=averageMatrix(simMatrix,1);
% syll_variability=mean(soundCorrVector)
% figure;plot(spikeCorrVector,soundCorrVector,'.');
% [R,P] = corrcoef(spikeCorrVector,soundCorrVector) 
