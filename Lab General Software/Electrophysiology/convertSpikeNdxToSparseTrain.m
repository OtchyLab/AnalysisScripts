function spikes = convertSpikeNdxToSparseTrain(spikeNdx, startNdx, endNdx)
%Creates a binary vector spike train as a sparse matrix.  

%spikeNdx are the indices of the spikes in arbitrary units.  If they are
%passed from motifData they represent the number of sample from the
%beginning of the motif at which the spike occured.

%The spike train will begin at startNdx and end at endNdx.  Any spikes
%outside of this window will not be represented in the train.

spikeNdx = spikeNdx(find(spikeNdx >= startNdx & spikeNdx <= endNdx));
spikeNdx = spikeNdx - startNdx + 1;
spikes = sparse(spikeNdx, ones(length(spikeNdx),1), ones(length(spikeNdx),1), endNdx-startNdx+1,1,length(spikeNdx));