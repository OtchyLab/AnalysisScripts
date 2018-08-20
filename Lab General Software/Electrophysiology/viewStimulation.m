function cell = viewStimulation(cell)
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
    zoom on;
    xlabel('Peri Stim Time (seconds)');
    ylabel('Amplified Voltage');
    title(['Cell ', num2str(cell.num), ' Stim File Num: ', num2str(stimFileNum)]); 

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
end

bQuest = input('Do you want to enter antidromic info? (1,0)');
if(bQuest)
	cell.antidromicIdentified = input('Was it identified?(y = 1 or n = 0)?');
	cell.antidromicAproxCurr = input('What stimulation current was used (uAmps)?');
	cell.antidromicGuessLat = input('What is the approximate latency (ms)?');
	cell.antidromicComments = input('Antidromic comments?','s');
	saveCell(cell);
end