function cell = viewStimForAntiFigure(cell)
sampRate = cell.exper.desiredInSampRate;
stimThreshold = 6.0; %Volts
preStim = .010; %seconds
postStim = .030; %seconds
preStimNdx = round(-preStim * sampRate);
postStimNdx = round(postStim * sampRate);

for(stimFileNum = cell.stimFileNums)
    %Find and extract stims
    sig = loadData(cell.exper, stimFileNum, cell.chanNum);
    aboveThres = find(sig > stimThreshold);
    downCrossNdx = find(diff(aboveThres) > sampRate/4);
    upCrossNdx = downCrossNdx + 1;
    downCross = [aboveThres(downCrossNdx); aboveThres(end)];
    upCross = [aboveThres(1); aboveThres(upCrossNdx)];
    stimNdx = repmat(upCross,1,postStimNdx-preStimNdx+1) + repmat([preStimNdx:postStimNdx],length(upCross),1);  
    %Becareful because last stim in file can get cut off.
    if(max(max(stimNdx)) > length(sig))
        stimNdx = stimNdx(1:end-1,:);
    end
    if(min(min(stimNdx)) < 1)
        stimNdx = stimNdx(2:end,:);
    end    
    stims = sig(stimNdx);
    
    h = figure(1005);
    plot([preStimNdx:postStimNdx]/sampRate,stims);         
    pause;
end
