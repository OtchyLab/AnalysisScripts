function spikes = convertSpikeTimesToTrain(spikeTime, startTime, endTime, binSize)
%Creates a binary vector spike train.
%spikeTime is the time of each spike in seconds.
%startTime is the time as with the vector begins, in seconds.  Spikes before 
%startTime are not used.
%endTime is the time as with the vector ends, in seconds.  Spikes after 
%endTime are not used.
%The size of the bins in seconds.

%NOTE:  The function will not work if you sampling at greater than 1MHz.

eps = 1/10000000;
spikeTime = spikeTime(find(spikeTime>=startTime & spikeTime<=endTime));
spikeNdx = floor((spikeTime - startTime + eps) / binSize + 1);
endNdx = floor((endTime - startTime + eps)/ binSize + 1);
spikes = zeros(endNdx, 1);
spikes(spikeNdx) = 1;


