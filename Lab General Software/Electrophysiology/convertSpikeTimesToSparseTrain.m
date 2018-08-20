function spikes = convertSpikeTimesToSparseTrain(spikeTime, startTime, endTime, binSize)
%Creates a binary vector spike train as a sparse matrix.  
%spikeTime is the time of each spike in seconds.
%startTime is the time as with the vector begins, in seconds.  Spikes before 
%startTime are not used.
%endTime is the time as with the vector ends, in seconds.  Spikes after 
%endTime are not used.
%The size of the bins in seconds.

