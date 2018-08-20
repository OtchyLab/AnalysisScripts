function alignmentPoints = alignTimes(audio_ts,video_ts) 
    numTimePoints = length(audio_ts); % don't need this? ideally we wont 
    % ever get to the end 
    numVidPoints = length(video_ts);
    curTimePoint = 1;
    alignmentPoints = zeros(1, numVidPoints);
    tol = 0.00009;
    for i = 1:numVidPoints
        if(audio_ts(numTimePoints) < video_ts(i))
            fprintf("Weird!");
            alignmentPoints(i) = NaN;
        else
            audVal = audio_ts(curTimePoint);
            vidVal = video_ts(i);
            while (abs(audVal-vidVal) > tol)
                curTimePoint = curTimePoint + 1; 
                audVal = audio_ts(curTimePoint);
            end
            timeVal = audio_ts(curTimePoint);
            alignmentPoints(i) = curTimePoint;
            curTimePoint = curTimePoint + 1; 
            while (i ~= numVidPoints && timeVal == audio_ts(curTimePoint))
                curTimePoint = curTimePoint + 1;
            end
        end
        
    end
end

