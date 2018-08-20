function [bSuccess, startNdx, endNdx, identity] = identifySyllablesInAudio(audio, sampRate, startThres, endThres, soundThres)

identity = [];

[startNdx, endNdx] = findSyllables(audio, sampRate, startThres, endThres, soundThres);

[bValid, startNdx, endNdx] = verifySyllablesSelection(audio, sampRate, startNdx, endNdx);

if(~bValid) 
    bSuccess = false;
    warning('Syllable identification failed.  Verify syllable function does not allow for correction.');
    return;
end


for(nSyll = 1:length(startNdx))
    %Display the sound...
    h = figure(1001);
    startWin = max(1, startNdx(nSyll) - (.05*sampRate));
    endWin = min(length(audio), endNdx(nSyll) + (.05*sampRate));    
    displayAudioSpecgram(audio(startWin:endWin), sampRate, startWin/sampRate, 8000, [-33,10], 10, .4);
    hold on;
    
    startNdx(nSyll) = max(1, startNdx(nSyll) - (.015*sampRate));
    endNdx(nSyll) = min(length(audio), endNdx(nSyll) + (.020*sampRate));    
    startTime = (startNdx(nSyll)-1) / sampRate;
    endTime = (endNdx(nSyll)-1) / sampRate;
    startLine = line([startTime, startTime], ylim, 'Color', 'yellow');
    endLine = line([endTime,endTime], ylim, 'Color', 'red');

    while(true)
        char = input('Instruction: 0-9, n,m,<,>,a', 's');
        if(char >= '0' & char <= '9')
            identity(nSyll) = str2num(char);
            break;
        elseif(char == 'd')
            identity(nSyll) = -1;
            break;
        elseif(char == 'n')
            delete(startLine);
            startNdx(nSyll) = max(1, startNdx(nSyll) - (.005*sampRate));
            startTime = (startNdx(nSyll)-1) / sampRate;
            startLine = line([startTime, startTime], ylim, 'Color', 'yellow');
        elseif(char == 'm')
            delete(startLine);
            startNdx(nSyll) = min(length(audio), startNdx(nSyll) + (.005*sampRate));
            startTime = (startNdx(nSyll)-1) / sampRate;
            startLine = line([startTime, startTime], ylim, 'Color', 'yellow');
        elseif(char == ',')
            delete(endLine);
            endNdx(nSyll) = max(1, endNdx(nSyll) - (.005*sampRate));
            endTime = (endNdx(nSyll)-1) / sampRate;
            endLine = line([endTime, endTime], ylim, 'Color', 'red');
        elseif(char == '.')
            delete(endLine);
            endNdx(nSyll) = min(length(audio), endNdx(nSyll) + (.005*sampRate));
            endTime = (endNdx(nSyll)-1) / sampRate;
            endLine = line([endTime, endTime], ylim, 'Color', 'yellow');
        end
    end
    close(h);
end
        
            
    
    
        





