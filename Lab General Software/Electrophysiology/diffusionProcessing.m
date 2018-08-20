function [stims, stimAvg, stimStd] = viewStimulation(exper, fileNumber, stimThreshold, chan)
sampRate = cell.exper.desiredInSampRate;
preStim = .010; %seconds
postStim = .050; %seconds
preStimNdx = round(-preStim * sampRate);
postStimNdx = round(postStim * sampRate);

if(~exist('chan'))
    chan = 1;
end

%Find and extract stims
sig = loadData(cell.exper, stimFileNum, chan);
aboveThres = find(sig < stimThreshold);
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
zoom on;
xlabel('Peri Stim Time (seconds)');
ylabel('Amplified Voltage');
title([' Stim File Num: ', num2str(stimFileNum)]); 

clear nStim;
while(true);
    pause(.05);
    char = get(h,'CurrentCharacter');
    set(h,'CurrentCharacter','~');
    if(char == ' ')
        break;     
    elseif(char == '.')
        if(~exist('nStim'))
            nStim = 1;
        elseif(nStim < size(stims,1))
            nStim = nStim + 1;
        end          
        v = axis; plot([preStimNdx:postStimNdx]/sampRate,stims(nStim,:));
        axis(v);
        zoom on;
        figure(h);
    elseif(char == ',')
        if(~exist('nStim'))
            nStim = 1;
        elseif(nStim >1)
            nStim = nStim - 1;
        end
        v = axis; plot([preStimNdx:postStimNdx]/sampRate,stims(nStim,:));
        axis(v);
        figure(h);
    elseif(char == 'a')
        v = axis; plot([preStimNdx:postStimNdx]/sampRate,stims);
        axis(v);
        figure(h);
    end
end        

avg = 