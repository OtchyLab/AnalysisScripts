function [R, startTime, endTime, motifMap] = cellMotifPairWiseCorrelations(cell, warpTemplate, bInstant, bRotate, gaussWidth, startMotif, endMotif, startTime,endTime)
%cell: the cell for which you want pairwise correlations.
%bWarp: should the spikes be time warped.
%bInstant: use instantaneous firing rate rather than spikes.
%gaussWidth: if non-zero, then filter signals with gaussian with 1/e point
%at gaussWidth seconds.
%R is the r values for each pair.  If a motif does not yet have a
%spikeTrain then the R value is not a NaN (not a number).  
%commonStartTime: the correlations are computed across the greatest common
%time interval for which all motif have spike data.  commonStartTime is the
%begining of this time interval.  If startTime and endTime are specifed the
%shorter signals are zero-padded.
%RValues, all values<1 in the R matrix, in the form of a vector.  If
%bRotate is true, then includes all values<1 in 100 R matrixies with each
%with different random rotations.

bWarp = length(warpTemplate) > 1;

if(~exist('startMotif'))
    startMotif =1;
    endMotif = length(cell.motifData);
end

sampRate = cell.exper.desiredInSampRate;

%Determine subset of times relative to the marker
%for which all motifs have spike data.
commonStartNdx = -20000000;
commonStopNdx = 20000000;
for(nMotif = startMotif:endMotif)
    motifData = cell.motifData{nMotif};
    if(isfield(motifData,'spikeNdx'))
        marker = motifData.marker;
        if(1-(marker-1) > commonStartNdx) disp(-1 * nMotif); end
        if((motifData.ndxStop-motifData.ndxStart+1) - (marker-1) < commonStopNdx) disp(nMotif); end
        commonStartNdx = max(commonStartNdx, 1-(marker-1));
        commonStopNdx = min(commonStopNdx, (motifData.ndxStop-motifData.ndxStart+1) - (marker-1));
    end
end
commonStartTime = ndx2time(commonStartNdx,sampRate);
commonEndTime = ndx2time(commonStopNdx,sampRate);

if(bWarp)
    startNdx = warpTemplate(1);
    endNdx = warpTemplate(end);
    startTime = ndx2time(startNdx,sampRate);
    endTime = ndx2time(endNdx,sampRate);
elseif(~exist('startTime'))
    startTime = commonStartTime;
    endTime = commonEndTime;
    startNdx = commonStartNdx
    endNdx = commonStopNdx;
else
    startNdx = startTime * sampRate + 1;
    endNdx = endTime * sampRate + 1;
end


%Remove motifs that have no spikes in time window.  %undefinied
%correlation.
nCount = 0;
motifMap = [];
for(nMotif = startMotif:endMotif)
    motifData = cell.motifData{nMotif};
    if(isfield(motifData,'spikeNdx') & length(find(motifData.spikeNdx-motifData.marker >= startNdx & motifData.spikeNdx-motifData.marker <=endNdx)))
        nCount = nCount + 1;
        motifMap(nCount) = nMotif;
    end
end


%Generate gauss for window
if(exist('gaussWidth') & gaussWidth > 0)
    sigma = gaussWidth / sqrt(2);
    x = [gaussWidth*-4:1/sampRate:gaussWidth*4];
    gauss = (1/sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));
end
    

%Clip out and filter the spike train of each motif
sigs = zeros(length(motifMap), endNdx - startNdx + 1);
for(nMotif = startMotif:endMotif)
    motifData = cell.motifData{nMotif};
    if(isfield(motifData,'spikeNdx') & length(find(motifData.spikeNdx-motifData.marker >= startNdx & motifData.spikeNdx-motifData.marker <=endNdx)))
        if(~bWarp)
            spikeTime = spikeTimeFromNdx(cell, nMotif);
        else
            spikeTime = getWarpedSpikeTimes(cell, nMotif, warpTemplate);
        end        
        
        if(bInstant)            
            %sig(nMotif,:) =  spikeTimes2InstantFiringRate(spike   
        else            
            sig = convertSpikeTimesToTrain(spikeTime, startTime, endTime, 1/sampRate);
        end
        
        if(bRotate)
            rotation = ceil(rand(1,1) * length(sig));
            sig = [sig(rotation:end);sig(1:rotation-1)];
        end        
        
        if(exist('gaussWidth') & gaussWidth > 0)
            sig  = conv(sig, gauss);
            sig = sig(floor(length(gauss)/2):end - floor((length(gauss)+1)/2));
        end
            
        sigs(find(motifMap == nMotif),:) = sig';
    end
end
clear motifData;

R = corrcoef(sigs');

%Compute correlations of all pairs for common subset of motif.
% sigsq = var(sigs');
% Q = zeros(length(cell.motifData)) * nan;
% biasedScale = 1/size(sigs,2);
% for(nMotif1 = 1:length(cell.motifData))
%     motifData1 = cell.motifData{nMotif1};
%     if(isfield(motifData1,'spikeNdx'))
%         for(nMotif2 = nMotif1:length(cell.motifData))
%             motifData2 = cell.motifData{nMotif2};
%             if(isfield(motifData2,'spikeNdx'))                
%                 Q(nMotif1,nMotif2) = cov(sigs(nMotif1,:),sigs(nMotif2,:))/
%             end
%         end
%     end
% end
        

%Perhaps faster?
%Given a matrix xcorr will complete all the pairwise correlations, but requires making        
%pairCorr = xcorr( allSpikeTrains,1,'unbiased');


