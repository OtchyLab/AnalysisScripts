function [maxCorr, motifShifts] = crossCorrelateBence(startMotif, endMotif, spikeTimes, sampRate, bInstant, bRotate, gaussWidth, startTime,endTime,sampRateOut)
         

%Aaron Andalman 2004
%Code for Bence base on cellMotifPairWiseCorrelations.
%The code is a bit difficult to read because it has remnants from the
%original version.

%INPUT
%spikeTimes: the cell array containing the spiketimes for each motif for which you want pairwise correlations.
%the spikestimes must be relative to some marker within the motif.
%sampRate: the sampling rate.
%bInstant: correlate instantaneous firing rate rather than spikes.
%gaussWidth: if non-zero, then filter signals with gaussian with 1/e point
%at gaussWidth seconds.
%startTime: spikes prior to this time are ignored.  Shorter trains are zero
%padded.  Choose this time based on your recordings and the length of the
%song.
%endTime: spikes after this time are ignored.
%sampRateOut: the sample rate of the fileted spike train
%reference: spike train or psth to which to correlate the spike times
%cellRef: the cell number that is to be used as reference
%OUTPUT
%motifShifts is an array containing the ideal shift in ms for each motif to best fit the reference spike train. 

%If bRotate is true, then includes all values<1 in 100 R matrixies with each
%with different random rotations.

% startMotif = 1;
% endMotif = length(spikeTimes);



startNdx = startTime * sampRate + 1;
endNdx = endTime * sampRate + 1;

%Remove motifs that have no spikes in time window.  
nCount = 0;
motifMap = [];
for(nMotif = startMotif:endMotif)
    %if(length(find(spikeTimes{nMotif} >= startTime & spikeTimes{nMotif} <= endTime)))
        nCount = nCount + 1;
        motifMap(nCount) = nMotif;
    %end
end

%Generate gauss for window
if(exist('gaussWidth') & gaussWidth > 0)
    sigma = gaussWidth / sqrt(2);
    x = [gaussWidth*-4:1/sampRate:gaussWidth*4];
    gauss = (1/sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));
end
    
%Clip out and filter the spike train of each motif
%sigs = zeros(length(motifMap), endNdx - startNdx + 1);%sometimes 2 instead of 1 here ???
maxCorr=ones(endMotif-startMotif+1);
maxCorr=-maxCorr; %set all correlation to -1 to start with;
for refMotif=startMotif:endMotif-1
    for sigMotif = refMotif+1:endMotif
        spikeTimeSig = spikeTimes{1,sigMotif};
        spikeTimeRef = spikeTimes{1,refMotif};
        %if(length(find(spikeTime >= startTime & spikeTime <= endTime)))      
            if(bInstant)       
                error('code not completed');
                %sig(nMotif,:) =  spikeTimes2InstantFiringRate(spike   
            else  
                ref = convertSpikeTimesToTrain(spikeTimeRef, startTime, endTime, 1/sampRate);
                sig = convertSpikeTimesToTrain(spikeTimeSig, startTime, endTime, 1/sampRate);
            end

            if(bRotate)
                rotationSig = ceil(rand(1,1) * length(sig));
                rotationRef = ceil(rand(1,1) * length(ref));
                sig = [sig(rotationSig:end);sig(1:rotationSig-1)];
                ref = [ref(rotationRef:end);sig(1:rotationRef-1)];
            end        

            if(exist('gaussWidth') & gaussWidth > 0)
                sig  = conv(sig, gauss);
                sig = sig(floor(length(gauss)/2):end - floor((length(gauss)+1)/2));
                sig = resample(sig,sampRateOut,sampRate)'; %resample to 500 Hz sampling
                ref  = conv(ref, gauss);
                ref = ref(floor(length(gauss)/2):end - floor((length(gauss)+1)/2));
                ref = resample(ref,sampRateOut,sampRate)'; %resample to 500 Hz sampling
            end
            [XCF, Lags, Bounds] = crosscorr(sig, ref, 5);
            [maxCorr(refMotif,sigMotif),I]=max(XCF);
            motifShifts(refMotif,sigMotif)=Lags(I)*(1000/sampRateOut); %shift in ms
    %         sigs(find(motifMap == nMotif),:) = sig';

               % sigs(nMotif-startMotif+1,:) = sig';
        %end
    end
 end



