function aud_ts = getAudTS(audData, frameRate)
    numPoints = length(audData);
    aud_ts = zeros(1, numPoints);
    time = 0; 
    increment = 1.0/frameRate; 
    for i = 1:numPoints
        aud_ts(i) = time;
        time = time + increment;
    end
end