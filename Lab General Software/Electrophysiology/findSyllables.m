function [syllableStart,syllableEnd] = findSyllables(audio, sampRate, startThres, endThres, soundThres)
%Returns start and end of syllables as array indices in audio file.

if(~exist('startThres'))
    startThres = -14;
end
if(~exist('endThres'))
    endThres = -15;
end

boxFiltLen = round(.015*sampRate);
boxFilter = repmat(1/boxFiltLen, boxFiltLen,1);

audio = audio - mean(audio);
audioLogPow = log(audio.^2 + eps);
audioLogPow = conv(audioLogPow, boxFilter);
audioLogPow = audioLogPow(length(boxFilter):end);

[upCross, junk] = detectThresholdCrossings(audioLogPow, startThres, true);
[junk, downCross] = detectThresholdCrossings(audioLogPow, endThres, true);

nSyll = 0;
lastEnd = 0;
while(true)
    ndxStart = find(upCross > lastEnd);
    if(length(ndxStart) == 0)
        break;
    end
    nSyll = nSyll + 1;
    syllableStart(nSyll) = upCross(ndxStart(1));
    ndxEnd = find(downCross > syllableStart(nSyll));
    if(length(ndxEnd) == 0)
        syllableEnd(nSyll) = length(audio); 
    end
    syllableEnd(nSyll) = downCross(ndxEnd(1));
    lastEnd = syllableEnd(nSyll);
end

for(nSyll = 1:length(syllableStart))
    ndx = find(audioLogPow(1:syllableStart(nSyll))< soundThres);
    if(length(ndx)>0)
        syllableStart(nSyll) = ndx(end);
    end
    ndx = find(audioLogPow(syllableEnd(nSyll):end) < soundThres);
    if(length(ndx)>0)
        syllableEnd(nSyll) = syllableEnd(nSyll) + ndx(1);
    end
end

%plot(audioLogPow);
%for(nSyll = 1:length(syllableStart))
%    line([syllableStart(nSyll), syllableStart(nSyll)], ylim, 'Color', 'red');
%    line([syllableEnd(nSyll), syllableEnd(nSyll)], ylim, 'Color', 'red');
%end
    



