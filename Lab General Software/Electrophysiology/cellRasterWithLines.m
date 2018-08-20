function numRows = cellRasterWithLines(cell, warpTemplate, startMotif, endMotif, highlight, startTime, endTime, color, startCount)
bWarp = length(warpTemplate) > 1;

sampRate = cell.exper.desiredInSampRate;
motifData = cell.motifData{1};
timeBuff = .1;
%startTime = - motifData.marker / sampRate;
%endTime = startTime + (motifData.ndxStop - motifData.ndxStart)/sampRate; 

%Plot raster of spikes from all motifs
hold on;
count = 0;
if(exist('startCount'))
    count = startCount;
end

if(~exist('color'))
    color = 'black';
end
   
commonStartNdx = -20000000;
commonStopNdx = 20000000;
for(nMotif = startMotif:endMotif)
    %Convert raw spike sample numbers in motifData.spikeNdx to
    %time in seconds relative to marker.
    motifData = cell.motifData{nMotif};
    if(isfield(motifData, 'spikeNdx'))
        count = count + 1;
        marker = motifData.marker;
        
        if(~bWarp)
            spikeTimeRelMarker = (motifData.spikeNdx - motifData.marker) / sampRate;
        else
            spikeTimeRelMarker = getWarpedSpikeTimes(cell, nMotif, warpTemplate);
        end
        
        for(nSpike = 1:length(spikeTimeRelMarker))
            l = line([spikeTimeRelMarker(nSpike),spikeTimeRelMarker(nSpike)], [count-.4,count+.4]);
            set(l,'Color',color);    
            if(exist('highlight') & ismember(nMotif,highlight))
                set(l,'Color','red');   
            end
                
            set(l,'LineWidth',1);
        end
        
        commonStartNdx = max(commonStartNdx, 1-(marker-1));
        commonStopNdx = min(commonStopNdx, (motifData.ndxStop-motifData.ndxStart+1) - (marker-1));
    end
end

if(bWarp)
    startTime = ndx2time(warpTemplate(1),sampRate);
    endTime = ndx2time(warpTemplate(end),sampRate);
elseif(~exist('startTime'))
    startTime = commonStartNdx/sampRate
    endTime = commonStopNdx/sampRate
end

xlim([startTime, endTime]);
xlabel('time secs');
ylabel('motif number');
ylim([0,count]);

numRows = count;