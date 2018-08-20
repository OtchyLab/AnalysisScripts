function cell = getAntidromicLatency(cell)
sampRate = cell.exper.desiredInSampRate;
stimThreshold = -5.0; %Volts
preStim = .010; %seconds
postStim = .030; %seconds
preStimNdx = round(-preStim * sampRate);
postStimNdx = round(postStim * sampRate);

latency = [];
notfired = 0;

h = figure(1005);

for(stimFileNum = cell.stimFileNums)
    %Find and extract stims
    sig = loadData(cell.exper, stimFileNum, cell.chanNum);
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
    
    figure(h);
    plot([preStimNdx:postStimNdx]/sampRate,stims);
    zoom on;
    xlabel('Peri Stim Time (seconds)');
    ylabel('Amplified Voltage');
    title(['Cell ', num2str(cell.num), ' Stim File Num: ', num2str(stimFileNum)]); 

    clear nStim;
    while(true);
        pause(.05);
        [x,y,char] = ginput(1);
        if(char == ' ')
            break;
        elseif(char == 'z')
            figure(h);
            zoom on;
            set(h,'CurrentCharacter','z');
            while(true)
                pause(.05);
                char = get(h,'CurrentCharacter');
                if(char == ' ')
                    break;
                end
            end
            figure(h);            
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
            plot([preStimNdx:postStimNdx]/sampRate,stims);
            figure(h);
        elseif(char == '2')
            beep; 
            notfired = notfired + 1;
        else
            latency = [latency, x];
        end
    end                 
end

cell.antiLatencies = latency;
cell.antiFailed = notfired;
saveCell(cell);