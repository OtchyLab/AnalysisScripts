function [rmsSyllab1, rmsSyllab3, rmsEnd] = computeMotifRMSJitter(cell)
    
sampRate = cell.exper.desiredInSampRate;

onset3 = [];
onset1 = [];
offset3 = [];
test = [];
for(nMotif = 1:length(cell.motifData))
    if(isfield(cell.motifData{nMotif},'spikeNdx'))
        motifData = cell.motifData{nMotif}
        warpMarks = motifData.warpMarks;
        marker = motifData.marker;
        if(marker ~= warpMarks(3))
            error('What?');
        end
        onset3 = [onset3, ndx2time(warpMarks(5)-(marker-1),sampRate)];
        onset1 = [onset1, ndx2time(warpMarks(1)-(marker-1),sampRate)];
        offset3 = [offset3, ndx2time(warpMarks(6)-(marker-1),sampRate)];
        test = [test, ndx2time(warpMarks(3)-(marker-1),sampRate)];
    end
end

rmsSyllab1 = std(onset1)
rmsSyllab3 = std(onset3)
rmsEnd = std(offset3)
std(test)