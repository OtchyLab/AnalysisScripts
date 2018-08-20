function checkCellMotifAlignment(cell, templateStartTime, templateEndTime)

if(isfield(cell,'motifData'))
    sampRate = cell.exper.desiredInSampRate;
    motifData = cell.motifData{1};
    
    %Plot the first motif as template of sorts in subplot 1.
    timeBuff = 0;
    h = figure(1003);
    subplot(5,1,1);
    startTime = - motifData.marker / sampRate;
    if(~exist('templateStartTime'))
        templateStartTime = startTime;
    end
    if(~exist('templateEndTime'))
        templateEndTime = startTime + (motifData.ndxStop - motifData.ndxStart)/sampRate; 
    end
    displayAudioSpecgram(motifData.audio, sampRate, startTime);
    xlim([templateStartTime-timeBuff, templateEndTime + timeBuff]);
    title(['cell:', num2str(cell.num), ' motif: 1']);
    
    %Show 4 motifs at a time in subplot 2 through 5.
    for(nMotif = 1:length(cell.motifData))
        subplot(5,1,mod(nMotif-1,4) + 2);
        motifData = cell.motifData{nMotif};
        startTime = -motifData.marker / sampRate;
        displayAudioSpecgram(motifData.audio, sampRate, startTime);
        xlim([templateStartTime-timeBuff, templateEndTime + timeBuff]);
        if(isfield(motifData,'spikeNdx'))
            title(['cell:', num2str(cell.num), ' motif:', num2str(nMotif)]);
        else
             title(['NOT INCLUDED cell:', num2str(cell.num), ' motif:', num2str(nMotif)]);
        end
        if(mod(nMotif,4) == 0)
            figure(1003);
            set(h,'CurrentCharacter','z');
            while(true)
                pause(.05);
                char = get(h,'CurrentCharacter');
                if(char == ' ')
                    break;
                end
            end
        end
    end
end